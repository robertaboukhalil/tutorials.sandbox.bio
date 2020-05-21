<script>
// TODO:
// - Tell user if answer is correct
// - Move to next lesson
// - Dropdown to browse other lessons
// - Initialize lessonNb state based on localStorage (i.e. where the user left off)
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
let bedtools = new Aioli("bedtools/2.29.2");

// Current lesson
let lesson = {};
let lessonNb = 0;

// Bed Files
let bedUser = { name: "Yours", contents: "" };
let bedUsage = { name: "Usage", contents: "Loading..." };

// UI State
let uiReady = false;
let uiError = "";
let uiInfo = "";
let uiCmd = "Loading...";


// -----------------------------------------------------------------------------
// Reactive statements
// -----------------------------------------------------------------------------

// Logic to update the lesson
$: {
	lesson = Lessons[lessonNb];
	lesson.goal.type = "goal";
	if(lessonNb > 0) {
		uiInfo = `<strong>Lesson Goal</strong>: Enter a <kbd>bedtools ${lesson.tool}</kbd> command ${lesson.description}:`;
		uiCmd = lesson.command;
		init(lesson.inputs);
	}
}

// Check whether user output is correct
$: bedUser.type = bedUser.contents == lesson.goal.contents ? "correct" : "incorrect";


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

/**
 * Initialize BED files for this lesson
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
	await bedtools.fs("chdir", files[0].directory);

	// Get documentation for relevant bedtools command
	bedUsage.contents = (await bedtools.exec(lesson.usage)).stderr;
}

/**
 * Run a bedtools command
 * @param {String} program Should always be "bedtools" 
 * @param {String} parameters Space-separated parameters to send to bedtools
 */
async function run(program, parameters)
{
	// Only accept bedtools commands
	uiError = "";
	if(program != "bedtools") {
		uiError = "Only bedtools commands are accepted.";
		return;
	}

	// Run bedtools with the parameters provided
	uiReady = false;
	let out = await bedtools.exec(parameters);

	// Save stdout/stderr
	bedUser.error = out.stderr != "";
	bedUser.contents = (out.stdout + out.stderr).trim();
	uiReady = true;
}


// -----------------------------------------------------------------------------
// On page load
// -----------------------------------------------------------------------------

onMount(async () => {
	// Initialize bedtools and output version
	await bedtools.init();
	await bedtools.exec("--version").then(d => console.log(d.stdout));

	uiReady = true;
	lessonNb = 1;
});


// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

<style>
.bg-light {
	background-color: rgb(232, 240, 240) !important;
}
</style>

<nav class="navbar navbar-expand-lg navbar-light bg-light">
	<div class="container">
		<span class="navbar-brand">&#x1F9EC; Bedtools Sandbox</span>
		<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNavDropdown" aria-controls="navbarNavDropdown" aria-expanded="false" aria-label="Toggle navigation">
			<span class="navbar-toggler-icon"></span>
		</button>
		<div class="collapse navbar-collapse border-left pl-2 pr-0" id="navbarNavDropdown">
			<ul class="navbar-nav">
				<li class="nav-item dropdown">
					<div class="dropdown-menu" aria-labelledby="navbarDropdownMenuLink">
						{#each Lessons as linkout, i}
							{#if i > 0}
								<a class="dropdown-item" href="#" on:click={() => lessonNb = i}><strong>Lesson {i}:</strong> {linkout.title}</a>
							{/if}
						{/each}
					</div>

					{#if lessonNb > 0}
					<a class="nav-link dropdown-toggle" href="#" id="navbarDropdownMenuLink" role="button" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
						<strong>Lesson {lessonNb}:</strong> {lesson.title}
					</a>
					{/if}
				</li>
			</ul>
		</div>
	</div>
</nav>


<main role="main">
	<div class="container mt-3">
		<!-- bedtools CLI -->
		<CommandLine
			info={uiInfo}
			error={uiError}
			command={uiCmd}
			on:execute={d => run(d.detail.program, d.detail.parameters)}
			disabled={!uiReady} />
 
		<!-- Visualize .bed files -->
		<div class="row mt-2">
			<div class="col-12">
				<BedViz beds={[ ...lesson.inputs, lesson.goal, bedUser]} />
			</div>
		</div>

		<!-- Inputs and outputs -->
		<div class="row mt-2">
			<div class="col-6">
				<Tabs tabs={[ bedUsage, ...lesson.inputs, lesson.goal ]} />
			</div>

			<div class="col-6">
				<Tabs tabs={[ bedUser ]} />
			</div>
		</div>
	</div>
</main>
