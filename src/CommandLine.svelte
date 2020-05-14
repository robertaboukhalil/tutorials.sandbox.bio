<script>
// Exports
export let command = "";        // Command to execute (e.g. samtools --version)
export let disabled = false;    // Whether to disable the input or not
export let preload = [];        // Which Aioli tools to preload

// Imports
import { onMount, createEventDispatcher } from "svelte";
import { Aioli } from "../public/aioli.js";  // FIXME: import { Aioli } from "@biowasm/aioli";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

const dispatch = createEventDispatcher();
// Aioli objects (key = tool name)
let Tools = {};


// -----------------------------------------------------------------------------
// Reactive Statements
// -----------------------------------------------------------------------------

// Split program name from parameters
$: program = command.split(" ").shift();
$: parameters = command.replace(`${program} `, "");


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Initialize the Aioli Worker with a tool of interest
async function initializeAioli(tool)
{
    let toolName, toolVersion;
    [toolName, toolVersion] = tool.split("/");
    console.log(`[CommandLine] Initializing ${toolName} v${toolVersion}`);

    Tools[toolName] = new Aioli(tool);
	await Tools[toolName].init();
    console.log(`[CommandLine] Done initializing ${toolName} v${toolVersion}`);
}

// Execute a command
function run()
{
    dispatch("execute", {
        program: program,
        parameters: parameters
    });
}

// Execute command when press Enter
function handleKeydown(event)
{
    if(event.key == "Enter")
        run();
}


// -----------------------------------------------------------------------------
// On load
// -----------------------------------------------------------------------------

// Pre-load Aioli modules of interest
onMount(async () => {
    for(let tool of preload)
        initializeAioli(tool);
});


// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

<div class="row">
    <div class="col-12">
        <div class="input-group">
            <div class="input-group-prepend">
                <span class="input-group-text">$</span>
            </div>
            <input
                type="text" id="command" class="form-control form-control-lg" style="font-family:monospace;"
                bind:value={command} on:keydown={handleKeydown}
                disabled={disabled} autofocus>
            <div class="input-group-append">
                <button class="btn btn-lg btn-primary" on:click={run}>
                    Run
                </button>
            </div>
        </div>
    </div>
</div>
