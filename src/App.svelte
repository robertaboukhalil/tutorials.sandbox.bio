<script>
// TODO:
// - Add ability to show a hint to the user
// - Press Enter to go to next lesson
// - Intro lesson that introduces bedtools and what it is
// - Initialize lessonNb state based on localStorage (i.e. where the user left off)
// - Test in Chrome
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
// TODO: load this from localStorage
let lessonHistory = {
	"lesson-id": {
		success: false,
		answer: "<user input>"
	}
};

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
	bedUser.contents = "";
	if(lessonNb > 0) {
		uiInfo = `<strong>Lesson Goal</strong>: Enter a <code>bedtools ${lesson.tool}</code> command ${lesson.description}:`;
		uiCmd = lesson.command;
		init(lesson.inputs);
	}
}

// TODO: save this to localStorage
$: console.log(lessonHistory);


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Initialize BED files for this lesson.
// Format: [{ name: "test.bed", contents: "chr1\t123\t456" }, ...]
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

	// If user was here before, show their answer
	// FIXME:
	// if(lessonHistory[lesson.id]) {
	// 	uiCmd = `bedtools ${lessonHistory[lesson.id].answer}`;
	// 	run("bedtools", lessonHistory[lesson.id].answer, true);
	// }

	// Get documentation for relevant bedtools command
	bedUsage.contents = (await bedtools.exec(lesson.usage)).stderr.trim();
}

// Run a bedtools command (input = CommandLine component message)
// Format: { program: "bedtools", args: "intersect", done: callback_fn() }
async function run(cli)
{
	// Only accept bedtools commands
	uiError = "";
	if(cli.program != "bedtools") {
		uiError = `Invalid command <kbd>${cli.program}</kbd>. Only <code>bedtools</code> commands are accepted.`;
		cli.done();
		return;
	}

	// Run bedtools with the parameters provided
	let out = await bedtools.exec(cli.args);

	// Save stdout/stderr
	bedUser.error = out.stderr != "";
	bedUser.contents = (out.stdout + out.stderr).trim();

	// Check whether user output is correct
	let success = bedUser.contents == lesson.goal.contents
	bedUser.type = success ? "correct" : "incorrect";
	lessonHistory[lesson.id] = {
		success: success,
		answer: cli.args
	};

	// If answer is correct, let the user know
	// FIXME: && !silent
	if(success) {
		jQuery("#modalSuccess").modal({ focus: false });  // don't focus on modal
	}

	cli.done();
}


// -----------------------------------------------------------------------------
// On load
// -----------------------------------------------------------------------------

onMount(async () => {
	// Initialize bedtools and output version
	await bedtools.init();
	await bedtools.exec("--version").then(d => console.log(d.stdout));

	// Prep UI
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

<!-- Header -->
<nav class="navbar navbar-expand-lg navbar-light bg-light pt-4 pb-4">
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
								<button
									class="btn btn-link dropdown-item pb-2 pt-2"
									style="vertical-align: baseline;"
									on:click={() => lessonNb = i}>
									<i class="fas fa-check" style="color: {lessonHistory[linkout.id] && lessonHistory[linkout.id].success ? "#3BA99C" : "#CCCCCC"}"></i>
									&nbsp;
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
		<!-- bedtools CLI -->
		<CommandLine
			info={uiInfo}
			error={uiError}
			command={uiCmd}
			disabled={!uiReady}
			on:execute={d => run(d.detail)} />
 
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
				<button type="button" class="btn btn-info" data-dismiss="modal" on:click={() => lessonNb++}>Go to the next lesson</button>
			</div>
		</div>
	</div>
</div>
