<script>
// Exports
export let command = "";        // Command to execute (e.g. samtools --version)
export let info = "";           // Info message to show above CLI
export let error = "";          // Error message to show below CLI
export let disabled = false;    // Whether to disable the input or not

// Imports
import { createEventDispatcher } from "svelte";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

let dispatch = createEventDispatcher();  // used to communicate with parent component
let program = null;                      // e.g. bedtools
let args = null;                         // e.g. intersect -a a.bed -b b.bed


// -----------------------------------------------------------------------------
// Reactive Statements
// -----------------------------------------------------------------------------

// Split program name from args
$: program = command.split(" ").shift();
$: args = command.replace(`${program} `, "");


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Execute a command
function run()
{
	// Send a message to parent component about what to execute
	disabled = true;
	dispatch("execute", {
		program: program,
		args: args,
		done: () => disabled = false
	});
}


// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

<style>
input {
	font-family: monospace;
}
</style>

<!-- Info message -->
<div class="row">
	<div class="col-12">
		<span class="text-muted">
			<span class="text-info">{@html info}&nbsp;</span>
		</span>
	</div>
</div>

<!-- CLI -->
<div class="row mt-2">
	<div class="col-12">
		<div class="input-group">
			<div class="input-group-prepend">
				<span class="input-group-text">$</span>
			</div>
			<input
				type="text"
				class="form-control form-control-lg"
				disabled={disabled}
				bind:value={command}
				on:keydown={event => event.key == "Enter" ? run() : null}
			/>
			<div class="input-group-append">
				<button
					class="btn btn-lg btn-info"
					disabled={disabled}
					on:click={run}
				>
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
