<script>
export let tabs = [];
export let active = tabs.length == 0 ? "" : tabs[0].name;
</script>

<style>
div {
    max-height: 350px;
    overflow: scroll;
}

.nav-item {
    font-size: 20px;
}
</style>

<ul class="nav nav-tabs mb-2">
    {#each tabs as tab}
        <li class="nav-item">
            <a class="nav-link" class:active={active == tab.name} href="#" on:click={() => active = tab.name}>{tab.name}</a>
        </li>
    {/each}
</ul>

<div>
{#each tabs as tab}
    {#if active == tab.name}
        {#if tab.error}
            <code>{tab.contents}</code>
        {:else}
            <pre>{tab.contents}</pre>
        {/if}
    {/if}
{/each}
</div>
