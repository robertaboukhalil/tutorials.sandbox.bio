<script>
import { Tools } from "../utils.js";
</script>

<style>
.bg-light {
	background-color: rgb(232, 240, 240) !important;
}
</style>

<nav class="navbar navbar-expand-lg navbar-light bg-light pt-4 pb-4">
    <div class="container">
        <span class="navbar-brand">&#x1F9EC; Sandbox.bio Tutorials</span>
        <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNavDropdown" aria-controls="navbarNavDropdown" aria-expanded="false" aria-label="Toggle navigation">
            <span class="navbar-toggler-icon"></span>
        </button>
    </div>
</nav>

<main role="main">
    <div class="container mt-3">
        <h4 class="mb-4">&nbsp;Tutorials:</h4>
        {#each Tools as tool}
        <div class="row mt-2 ml-1">
            <a
                class="btn btn-outline-primary pb-2 pt-2"
                href="?tool={tool.id}">
                <strong>{tool.id}:</strong> {tool.description}
            </a>
        </div>
        {/each}
    </div>
</main>
