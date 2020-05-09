<script>
import { onMount } from "svelte";
import { Aioli } from "../public/aioli.js";
// FIXME: Import from npm once push Aioli fixes
// import { Aioli } from "@biowasm/aioli";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

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
let Output = "Loading...";


// -----------------------------------------------------------------------------
// Reactive Statements
// -----------------------------------------------------------------------------

// $: Bedtools.exec("intersect -a ")

$: bedtoolsIntersect(Bed1, Bed2);

async function bedtoolsIntersect(bed1, bed2)
{
	// Don't do anything until Aioli is done initializing
	if(!Bedtools.ready)
		return;

	// Generate Blob objects from strings
	let blob1 = new Blob([ bed1 ], { type: "text/plain" }),
		blob2 = new Blob([ bed2 ], { type: "text/plain" });

	// Mount those Blobs as files on the virtual file system
	let f1 = await Aioli.mount(blob1, "a.bed");
	let f2 = await Aioli.mount(blob2, "b.bed");

	// Run bedtools insersect
	let cmd = `intersect -a ${f1.path} -b ${f2.path}`;
	console.log(cmd);
	let output = await Bedtools.exec(cmd);
	console.warn(output.stderr);
	console.log(output.stdout);

	Output = output.stdout;
}


// -----------------------------------------------------------------------------
// On page load
// -----------------------------------------------------------------------------

// Initialize bedtools2 and output version
onMount(async () => {
	Bedtools.init()
		.then(() => Bedtools.exec("--version"))
		.then(d => {
			console.warn(`Loaded ${d.stdout}`)
			bedtoolsIntersect(Bed1, Bed2);
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
		<div class="row">
			<!-- Input -->
			<div class="col-md-4 mb-2">
				<h4 class="mb-3">Bed File A</h4>
				<textarea class="form-control" rows="10" bind:value={Bed1}></textarea>
			</div>
			<div class="col-md-4 mb-2">
				<h4 class="mb-3">Bed File B</h4>
				<textarea class="form-control" rows="10" bind:value={Bed2}></textarea>
			</div>

			<!-- Result -->
			<div class="col-md-4 mb-2">
				<h4 class="mb-3">Result</h4>
				<textarea class="form-control" rows="10" bind:value={Output}></textarea>
			</div>
		</div>
	</div>
</main>
