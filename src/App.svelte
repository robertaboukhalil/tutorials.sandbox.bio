<script>
import { onMount } from "svelte";
import { Aioli } from "@biowasm/aioli";

import { Lessons } from "./lessons.js";
import BedViz from "./BedViz.svelte";
import CommandLine from "./CommandLine.svelte";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

let Filename1 = "a.bed";
let Filename2 = "b.bed"

let Bedtools = new Aioli("bedtools/2.29.2");
let Bed1 = `
chr1	10	20	a1	1	+
chr1	100	200	a2	2	-
`.trim();
let Bed2 = `
chr1	20	30	b1	1	+
chr1	90	101	b2	2	-
chr1	100	110	b3	3	+
chr1	200	210	b4	4	+
`.trim();
let StdOut = "";
let StdErr = "";
let Cmd = `bedtools intersect -a ${Filename1} -b ${Filename2}`;

let UI = {
	error: "",
	info: "Enter the <i>bedtools</i> command that will generate the desired output:"
};


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
	let directory = files[0].directory;  // this should always be /data, but just in case
	await Bedtools.fs("chdir", directory);
}

/**
 * Run a bedtools command
 * 
 * @param {String} program Should always be "bedtools" 
 * @param {String} parameters Space-separated parameters to send to bedtools
 * @return {Object} Result of the command: { stdout: "<stdout>", stderr: "<stderr>" }  # FIXME:
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
	let output = await Bedtools.exec(parameters);
	StdOut = output.stdout;
	StdErr = output.stderr.replace(/\n/g, '<br>');
	console.log(output);
}


// -----------------------------------------------------------------------------
// On page load
// -----------------------------------------------------------------------------

onMount(async () => {
	// Initialize bedtools and output version
	await Bedtools.init();
	await Bedtools.exec("--version").then(d => console.log(d.stdout));

	// Initialize BED files
	init([
		{ name: "a.bed", contents: Bed1 },
		{ name: "b.bed", contents: Bed2 },
	]);
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
			<h2 class="display-5">Lesson 1</h2>
			<p class="lead">Using bedtools intersect to do xyz</p>
		</div>
	</div>

	<div class="container">
		<!-- bedtools CLI -->
        <CommandLine
            command={Cmd}
            info={UI.info}
            error={UI.error}
            on:execute={d => run(d.detail.program, d.detail.parameters)}
        />
 
		<!-- Visualize .bed files -->
		<div class="row">
			<div class="col-12">
				<BedViz beds={[Bed1, Bed2, StdOut]} names={[Filename1, Filename2, "Your output"]} />
			</div>
		</div>

		<!-- TODO: Replace with Tab component -->
		<div class="row">
			<!-- Input -->
			<div class="col-md-4 mb-2">
				<h4 class="mb-3">a.bed</h4>
				<textarea class="form-control" rows="10" bind:value={Bed1}></textarea>
			</div>
			<div class="col-md-4 mb-2">
				<h4 class="mb-3">b.bed</h4>
				<textarea class="form-control" rows="10" bind:value={Bed2}></textarea>
			</div>

			<!-- Result -->
			<div class="col-md-4 mb-2">
				<h4 class="mb-3">Result</h4>
				<div id="output">
					{@html StdOut || "Loading ..."}
					<strong><span style="color: red;">
						{@html StdErr}
					</span></strong>
				</div>
			</div>
		</div>
	</div>
</main>
