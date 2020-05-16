<script>
// TODO:
// - Optimize reactive statements (some get triggered very often)
// - Tell user if answer is correct
// - Move to next lesson
// - Initialize UI.lesson state based on localStorage (i.e. where the user left off)
// - How does it look on mobile?

import { onMount } from "svelte";
import { Aioli } from "@biowasm/aioli";
import { Lessons } from "./lessons.js";

import Tabs from "./Tabs.svelte";
import BedViz from "./BedViz.svelte";
import CommandLine from "./CommandLine.svelte";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

// Bedtools WebAssembly module
let Bedtools = new Aioli("bedtools/2.29.2");

// Lesson in progress
let Lesson = {};
let BedUser = {
	name: "Yours",
	contents: "",
	color: "red",
	error: false
}

// Tabs to show user
let BedTabs = {
	input: {},
	output: {},
	usage: {
		name: "Usage",
		contents: "Loading..."
	},
};

// UI State
let UI = {
	lesson: 0,
	ready: false,
	error: "",
	info: "Enter the <i>bedtools</i> command that will generate the desired output:",
	cmd: "bedtools intersect -a a.bed -b b.bed"
};


// -----------------------------------------------------------------------------
// Reactive statements
// -----------------------------------------------------------------------------

// Logic to update the lesson
$: Lesson = Lessons[UI.lesson];
$: init(Lesson.inputs);
$: BedUser.color = BedUser.contents == Lesson.goal.contents ? "green" : "red";

// Update tabs to display
$: BedTabs.input = [
	BedTabs.usage,
	Lesson.goal,
	...Lesson.inputs
];

$: BedTabs.output = [BedUser];


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

	// Get usage of relevant bedtools command
	let usage = await Bedtools.exec(Lesson.usage);
	BedTabs.usage.contents = usage.stderr;
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
	if(program != "bedtools") {
		UI.error = "Only bedtools commands are accepted.";
		return;
	}

	// Run bedtools with the parameters provided
	UI.ready = false;
	let output = await Bedtools.exec(parameters);

	// Save stdout/stderr
	BedUser.error = output.stderr != "";
	BedUser.contents = BedUser.error ? output.stderr.trim() : output.stdout.trim();
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
			info={UI.info}
			error={UI.error}
			bind:command={UI.cmd}
			on:execute={d => run(d.detail.program, d.detail.parameters)}
			disabled={!UI.ready} />
 
		<!-- Visualize .bed files -->
		<div class="row">
			<div class="col-12">
				<BedViz beds={[ ...Lesson.inputs, Lesson.goal, BedUser ]} />
			</div>
		</div>

		<!-- Inputs and outputs -->
		<div class="row">
			<div class="col-6">
				<Tabs tabs={BedTabs.input} />
			</div>

			<div class="col-6">
				<Tabs tabs={BedTabs.output} />
			</div>
		</div>
	</div>
</main>
