<script>
// Exports
export let command = "";        // Command to execute (e.g. samtools --version)
export let disabled = false;    // Whether to disable the input or not
export let info = "";           // Info message to show above CLI
export let error = "";          // Error message to show below CLI

// Imports
import { createEventDispatcher } from "svelte";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

let dispatch = createEventDispatcher();  // used to communicate with parent component
let program = null;                      // e.g. bedtools
let parameters = null;                   // e.g. intersect -a a.bed -b b.bed


// -----------------------------------------------------------------------------
// Reactive Statements
// -----------------------------------------------------------------------------

// Split program name from parameters
$: program = command.split(" ").shift();
$: parameters = command.replace(`${program} `, "");


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Execute a command
function run()
{
	// Send a message to parent component about what to execute
	dispatch("execute", {
		program: program,
		parameters: parameters
	});
}


// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

<!-- Info message -->
<div class="row">
	<div class="col-12">
		<small class="text-muted">
			<span class="text-info">{@html info}&nbsp;</span>
		</small>
	</div>
</div>

<!-- CLI -->
<div class="row">
	<div class="col-12">
		<div class="input-group">
			<div class="input-group-prepend">
				<span class="input-group-text">$</span>
			</div>
			<input
				type="text"
				class="form-control form-control-lg"
				style="font-family: monospace"
				bind:value={command}
				on:keydown={event => event.key == "Enter" ? run() : null}
				disabled={disabled}
				autofocus
			>
			<div class="input-group-append">
				<button class="btn btn-lg btn-info" on:click={run}>
					Run
				</button>
			</div>
		</div>
	</div>
</div>

<!-- Error message -->
<div class="row">
	<div class="col-12">
		<small class="text-muted">
			<span class="text-danger">{@html error}&nbsp;</span>
		</small>
	</div>
</div>
