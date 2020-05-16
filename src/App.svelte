<script>
// TODO: 
// - Check if output matches expected answer
// - Initialize UI.lesson state based on localStorage (i.e. where the user left off)
// - Replace hardcoded textareas with Tab component

import { onMount } from "svelte";
import { Aioli } from "@biowasm/aioli";

import { Lessons } from "./lessons.js";
import BedViz from "./BedViz.svelte";
import CommandLine from "./CommandLine.svelte";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

// Bedtools WebAssembly module
let Bedtools = new Aioli("bedtools/2.29.2");
let StdOut = "Loading...";
let StdErr = "";

// Lesson in progress
let Lesson = {};
let BedUser = {
	name: "yours.bed",
	contents: ""
}

// UI State
let UI = {
	ready: false,
	lesson: 0,
	error: "",
	info: "Enter the <i>bedtools</i> command that will generate the desired output:"
};


// -----------------------------------------------------------------------------
// Reactive statements
// -----------------------------------------------------------------------------

// Logic to update the lesson
$: Lesson = Lessons[UI.lesson];
$: init(Lesson.inputs);


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

/**
 * Initialize BED files for this lesson
 * 
 * @param {Array} bedFiles List of BED files: [{ name: "test.bed", contents: "chr1\t123\t456" }, ...]
 * @return null
 */
async function init(bedFiles)
{
	// Generate Blob objects from strings and mount them as files
	let files = [];
	for(let fileInfo of bedFiles) {
		let blob = new Blob([ fileInfo.contents ], { type: "text/plain" });
		files.push(await Aioli.mount(blob, fileInfo.name));
	}

	// At this point, the files are mounted to /data. Let's update the working
	// directory so users can refer to the .bed files without specifying /data.
	if(files.length == 0)
		return;
	let directory = files[0].directory;
	await Bedtools.fs("chdir", directory);
}

/**
 * Run a bedtools command
 * 
 * @param {String} program Should always be "bedtools" 
 * @param {String} parameters Space-separated parameters to send to bedtools
 */
async function run(program, parameters)
{
	// Only accept bedtools commands
	UI.error = "";
	UI.ready = false;
	if(program != "bedtools") {
		UI.error = "Only bedtools commands are accepted.";
		return;
	}

	// Run bedtools with the parameters provided
	let output = await Bedtools.exec(parameters);
	BedUser.contents = output.stdout;

	// Display stdout/stderr
	StdOut = output.stdout;
	StdErr = output.stderr.replace(/\n/g, "<br>");

	UI.ready = true;
}


// -----------------------------------------------------------------------------
// On page load
// -----------------------------------------------------------------------------

onMount(async () => {
	// Initialize bedtools and output version
	await Bedtools.init();
	await Bedtools.exec("--version").then(d => console.log(d.stdout));

	UI.ready = true;
	UI.lesson = 1;
});


// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

<style>
textarea {
	font-family: Consolas, monospace;
	white-space: pre;
	word-wrap: normal;
}
#output {
	border: 1px solid lightgray;
	font-family: Consolas, monospace;
	padding: 2px;
	word-wrap: normal;
	height: 400px;
}
</style>

<nav class="navbar navbar-expand-md navbar-dark fixed-top bg-dark">
	<a class="navbar-brand" href="/">Bedtools Sandbox</a>
	<div class="collapse navbar-collapse">
		<ul class="navbar-nav mr-auto"></ul>
		<ul class="navbar-nav">
			<li class="nav-item active">
				<a class="nav-link" href="/">Code</a>
			</li>
		</ul>
	</div>
</nav>

<main role="main">
	<div class="jumbotron mt-4 pb-3">
		<div class="container">
			<h2 class="display-5">{Lesson.title}</h2>
			<p class="lead">{Lesson.description}</p>
		</div>
	</div>

	<div class="container">
		<!-- bedtools CLI -->
		<CommandLine
			command={Lesson.command}
			info={UI.info}
			error={UI.error}
			on:execute={d => run(d.detail.program, d.detail.parameters)}
			disabled={!UI.ready}
		/>
 
		<!-- Visualize .bed files -->
		<div class="row">
			<div class="col-12">
				<BedViz beds={[ ...Lesson.inputs, Lesson.goal, BedUser ]} />
			</div>
		</div>

		<!-- Inputs and outputs -->
		<div class="row">
			{#each Lesson.inputs as bedFile}
				<div class="col-md-4 mb-2">
					<h4 class="mb-3">{bedFile.name}</h4>
					<textarea class="form-control" rows="10" bind:value={bedFile.contents}></textarea>
				</div>
			{/each}

			<div class="col-md-4 mb-2">
				<h4 class="mb-3">Result</h4>
				<div id="output">
					{@html StdOut}
					<strong><span style="color: red;">
						{@html StdErr}
					</span></strong>
				</div>
			</div>
		</div>
	</div>
</main>
