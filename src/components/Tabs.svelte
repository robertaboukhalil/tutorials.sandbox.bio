<script>
// Exports
export let tabs = [];
export let active = tabs.length == 0 ? "" : tabs[0].name;
export let scroll = false;

// Imports
import { afterUpdate } from "svelte";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

let tab = null;
let contents = "";

// Once DOM is updated, scroll to the bottom of the div
afterUpdate(() => {
	if(scroll)
		tab.scrollTop = tab.scrollHeight;
});


// -----------------------------------------------------------------------------
// Reactive Statements
// -----------------------------------------------------------------------------

$: tabActive = tabs.filter(t => t.name == active).pop();
$: errors = tabActive.contents.match(/(.*error:.*)/gi);


// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

<style>
div {
	max-height: 350px;
	overflow: scroll;
}
</style>

<ul class="nav nav-tabs mb-0">
	{#each tabs as tab}
		<li class="nav-item">
				<a class="nav-link" href="#" class:active={active == tab.name} on:click={() => active = tab.name}>{tab.name}</a>
		</li>
	{/each}
</ul>

<div bind:this={tab} class="border-bottom border-right border-left border-default p-3">
	<pre>{tabActive.contents}<br /><br /></pre>
	{#if errors != "" && errors != null}
		<pre style="padding:10px; border:3px dotted red;">{@html errors.join("\n")}<br /><br /></pre>
	{/if}
</div>
