<script>
export let lessons = [];
export let toolID = "";
export let toolName = "";
export let toolAioli = "";
export let toolCLI = "";

import { onMount } from "svelte";
import { Aioli } from "@biowasm/aioli";
import { BedFile } from "../utils.js";

import Tabs from "./Tabs.svelte";
import VizBed from "./VizBed.svelte";
import CommandLine from "./CommandLine.svelte";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

// WebAssembly module
let aioli = new Aioli(toolAioli);

// Current lesson
let lesson = {};
let lessonNb = 0;
let lessonAnswers = {};

// Bed Files
let bedUser = new BedFile("Yours");
let bedUsage = new BedFile("Usage");
let bedGoal = new BedFile("Goal", "goal");

// UI State
let uiReady = false;
let uiError = "";
let uiInfo = "";
let uiCmd = "Loading...";
let uiBtnSuccess = null;


// -----------------------------------------------------------------------------
// Reactive statements
// -----------------------------------------------------------------------------

// Logic to update the lesson
$: lesson = lessons[lessonNb];
$: lessonNb > 0 ? init(lesson) : null;
$: lesson.answer = (lessonAnswers[lesson.id] || {}).answer;
$: lesson.answer = lesson.answer != lesson.tool ? lesson.answer : null;

// Get/set user's answers to localStorage to keep state when revisit the site
$: lessonAnswers = JSON.parse(localStorage.getItem(toolID) || "{}");
$: localStorage.setItem(toolID, JSON.stringify(lessonAnswers));

$: bedUser.type = bedUser.contents == bedGoal.contents ? "correct" : "incorrect";


// -----------------------------------------------------------------------------
// Aioli functions
// -----------------------------------------------------------------------------

// Execute a command
async function exec(args)
{
	if(!aioli.ready || args == null || args == "")
		return "";
	let result = await aioli.exec(args);
	return (result.stdout + result.stderr).trim();
}

// Mount an array of BedFile objects
async function mount(bedFiles)
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
	await aioli.fs("chdir", files[0].directory);
}


// -----------------------------------------------------------------------------
// Lesson functions
// -----------------------------------------------------------------------------

// Initialize a lesson
async function init(lesson)
{
	// Reset UI
	uiInfo = `<strong>Lesson Goal</strong>: Enter a <code>${toolCLI} ${lesson.tool}</code> command ${lesson.description}:`;
	uiCmd = `${toolCLI} ${lesson.answer || lesson.tool}`;
	bedUser.contents = "";

	// Mount the files we need in this lesson
	await mount(lesson.inputs);

    // Run tool
	[ bedGoal.contents, bedUser.contents, bedUsage.contents ] = await Promise.all([
		exec(lesson.goal),     // Run the command that should give the right output
		exec(lesson.answer),   // If user had previously entered an answer, run it
		exec(lesson.usage),    // Get relevant CLI documentation
	]);
}

// Execute a user's command obtained from CommandLine component
// Format: { program: "bedtools", args: "intersect", done: callback_fn() }
async function run(cli)
{
	// Only support commands for current tool
	uiError = "";
	if(cli.program != toolCLI) {
		uiError = `Invalid command <kbd>${cli.program}</kbd>. Please enter a <code>${toolCLI}</code> command.`;
		cli.done();
		return;
	}

	// Run the tool
	bedUser.contents = await exec(cli.args);
	let success = bedUser.contents == bedGoal.contents;

	// Save user input
	lessonAnswers[lesson.id] = {
		success: success,
		answer: cli.args
	};

	// If answer is correct, let the user know
	if(success) {
		jQuery("#modalSuccess").modal({ focus: false });
		setTimeout(() => uiBtnSuccess.focus(), 150);
	}
	cli.done();
}


// -----------------------------------------------------------------------------
// On load
// -----------------------------------------------------------------------------

onMount(async () => {
	// Initialize tool and output version
	await aioli.init();
	await aioli.exec("--version").then(d => console.log(d.stdout));

	// Prep UI
	uiReady = true;

	// Launch first unfinished lesson
	let i;
	for(i = 1; i < lessons.length; i++)
		if(!((lessonAnswers[lessons[i].id] || {}).success))
			break;
	lessonNb = i;
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

<!-- Header -->
<nav class="navbar navbar-expand-lg navbar-light bg-light pt-4 pb-4">
	<div class="container">
		<span class="navbar-brand">&#x1F9EC; {toolName} Sandbox</span>
		<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNavDropdown" aria-controls="navbarNavDropdown" aria-expanded="false" aria-label="Toggle navigation">
			<span class="navbar-toggler-icon"></span>
		</button>
		<div class="collapse navbar-collapse border-left pl-2 pr-0" id="navbarNavDropdown">
			<ul class="navbar-nav">
				<li class="nav-item dropdown">
					<div class="dropdown-menu" aria-labelledby="navbarDropdownMenuLink">
						{#each lessons as linkout, i}
							{#if i > 0}
								<button
									class="btn btn-link dropdown-item pb-2 pt-2 {linkout.id == lesson.id ? "bg-light" : ""}"
									style="vertical-align: baseline;"
									on:click={() => lessonNb = i}>
									<i class="fas fa-check" style="color: {(lessonAnswers[linkout.id] || {}).success ? "#3BA99C" : "#CCCCCC"}"></i>&nbsp;
									<strong>Lesson {i}:</strong> {linkout.title}
								</button>
							{/if}
						{/each}
					</div>

					{#if lessonNb > 0}
					<btn class="btn btn-link nav-link dropdown-toggle" id="navbarDropdownMenuLink" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
						<strong>Lesson {lessonNb}:</strong> {lesson.title}
					</btn>
					{/if}
				</li>
			</ul>
		</div>
	</div>
</nav>

<!-- Tutorial -->
<main role="main">
	<div class="container mt-3">
		<!-- CLI -->
		<CommandLine
			info={uiInfo}
			error={uiError}
			command={uiCmd}
			disabled={!uiReady}
			on:execute={d => run(d.detail)} />
 
		<!-- Visualize .bed files -->
		<div class="row mt-2">
			<div class="col-12">
				<VizBed beds={[ ...lesson.inputs, bedGoal, bedUser]} />
			</div>
		</div>

		<!-- Inputs and outputs -->
		<div class="row mt-2">
			<div class="col-6">
				<Tabs tabs={[ bedUsage, ...lesson.inputs, bedGoal ]} />
			</div>

			<div class="col-6">
				<Tabs tabs={[ bedUser ]} scroll="true" />
			</div>
		</div>
	</div>
</main>

<!-- Success Modal -->
<div id="modalSuccess" class="modal fade" tabindex="-1" role="dialog">
	<div class="modal-dialog" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<h5 class="modal-title">That's correct!</h5>
			</div>
			<div class="modal-footer">
				<button bind:this={uiBtnSuccess} type="button" class="btn btn-primary" data-dismiss="modal" on:click={() => lessonNb++}>Go to the next lesson</button>
			</div>
		</div>
	</div>
</div>
