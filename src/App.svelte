<script>
import { onMount } from "svelte";
import { Aioli } from "../public/aioli.js";
import BedViz from "./BedViz.svelte";

// FIXME: Import from npm once push Aioli fixes
// import { Aioli } from "@biowasm/aioli";


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


// -----------------------------------------------------------------------------
// Reactive Statements
// -----------------------------------------------------------------------------

// $: Bedtools.exec("intersect -a ")

$: bedtoolsIntersect(Cmd, Bed1, Bed2);

async function bedtoolsIntersect(cmd, bed1, bed2)
{
	// Don't do anything until Aioli is done initializing
	if(!Bedtools.ready)
		return;

	// Generate Blob objects from strings
	let blob1 = new Blob([ bed1 ], { type: "text/plain" });
	let blob2 = new Blob([ bed2 ], { type: "text/plain" });

	// Mount those Blobs as files on the virtual file system
	let f1 = await Aioli.mount(blob1, Filename1);
	let f2 = await Aioli.mount(blob2, Filename2);

	// Run bedtools insersect
	let final_command = cmd.replace("bedtools ", "")
													.replace(Filename1, f1.path)
													.replace(Filename2, f2.path)
													.trim();

	console.log("command:", final_command);
	let output = await Bedtools.exec(final_command);
	console.warn("stderr:", output.stderr);
	console.log("stdout:", output.stdout);
	StdOut = output.stdout;
	StdErr = output.stderr.replace(/\n/g, '<br>');
}


// -----------------------------------------------------------------------------
// On page load
// -----------------------------------------------------------------------------

// Initialize bedtools2 and output version
onMount(async () => {
	Bedtools.init()
		.then(() => Bedtools.exec("--version"))
		.then(d => {
			console.log(`Loaded ${d.stdout}`)
			bedtoolsIntersect(Cmd, Bed1, Bed2);
		});
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
			<h2 class="display-5">Title</h2>
			<p class="lead">Description</p>
		</div>
	</div>

	<div class="container">
			<!-- bedtools command -->
			<input type="text" class="form-control" bind:value={Cmd} />
			<BedViz beds={[Bed1, Bed2, StdOut]} names={[Filename1, Filename2, "Your output"]} />
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
