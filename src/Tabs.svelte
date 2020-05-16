<script>
export let tabs = [];
export let active = tabs.length == 0 ? "" : tabs[0].name;
</script>

<ul class="nav nav-tabs mb-3">
    {#each tabs as tab}
        <li class="nav-item">
            <a class="nav-link" class:active={active == tab.name} href="#" on:click={() => active = tab.name}>{tab.name}</a>
        </li>
    {/each}
</ul>

{#each tabs as tab}
    {#if active == tab.name}
        <pre>
            {tab.contents}
        </pre>
    {/if}
{/each}
