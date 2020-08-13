<script>
export let lessons = [];
export let toolID = "";
export let toolName = "";
export let toolAioli = "";
export let toolCLI = "";
export let postprocess = d => d;


// -----------------------------------------------------------------------------
// Imports
// -----------------------------------------------------------------------------

import { onMount } from "svelte";
import jQuery from "jquery";
import { Aioli } from "@biowasm/aioli";
import { SomeFile } from "../utils.js";

import Tabs from "./Tabs.svelte";
import CommandLine from "./CommandLine.svelte";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

// WebAssembly module
let aioli = new Aioli(toolAioli);
let tools = [
	{ id: "bedtools", description: "Wrangle .bed genomic interval files" },
	{ id: "bowtie", description: "Aligns reads to the genome" },
];

// Current lesson
let lesson = {};
let lessonNb = 0;
let lessonAnswers = {};

// Files displayed to the user
let fileUser = new SomeFile("Yours");
let fileUsage = new SomeFile("Usage");
let fileGoal = new SomeFile("Goal", "goal");

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

// Set the type to "correct" if user got it right (this is used by the viz)
$: fileUser.type = postprocess(fileUser.contents) == postprocess(fileGoal.contents) ? "correct" : "incorrect";


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

// Mount an array of SomeFile objects
async function mount(someFiles)
{
	// Generate Blob objects from strings and mount them as files
	let files = [];
	for(let fileInfo of someFiles) {
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
	fileUser.contents = "";

	// Mount the files we need in this lesson
	await mount(lesson.inputs);

    // Run tool
	[ fileGoal.contents, fileUser.contents, fileUsage.contents ] = await Promise.all([
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
	fileUser.contents = await exec(cli.args);
	let success = postprocess(fileUser.contents) == postprocess(fileGoal.contents);

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
			<ul class="navbar-nav w-100">
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
				<li class="nav-item dropdown ml-auto">
					<div class="dropdown-menu dropdown-menu-right" aria-labelledby="navbarDropdownMenuLink2">
						{#each tools as tool}
							<a
								class="btn btn-link dropdown-item pb-2 pt-2 {tool.id == toolID ? "bg-light" : ""}"
								style="vertical-align: baseline;"
								href="?tool={tool.id}">
								<strong>{tool.id}:</strong> {tool.description}
							</a>
						{/each}
					</div>

					<btn class="btn btn-link nav-link dropdown-toggle" id="navbarDropdownMenuLink2" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
						<strong>{toolID}</strong>
					</btn>
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
 
        <slot name="viz" files={[ ...lesson.inputs, fileGoal, fileUser]} />

		<!-- Inputs and outputs -->
		<div class="row mt-2">
			<div class="col-6">
				<Tabs tabs={[ fileUsage, ...lesson.inputs, fileGoal ]} />
			</div>

			<div class="col-6">
				<Tabs tabs={[ fileUser ]} scroll="true" />
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
