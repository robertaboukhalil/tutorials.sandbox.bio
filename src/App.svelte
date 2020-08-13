<script>
// Bootstrap
import popper from "popper.js";
import Bootstrap from "bootstrap";
import "bootstrap/dist/css/bootstrap.min.css";

// Home page
import Home from "./components/Home.svelte";

// CLI Tools
import Bedtools from "./apps/Bedtools.svelte";
import Bowtie from "./apps/Bowtie.svelte";

// Reactive Statements
$: tool = new URL(window.location).searchParams.get("tool");
</script>

{#if tool == "bedtools"}
	<Bedtools />
{:else if tool == "bowtie"}
	<Bowtie />
{:else}
	<Home />
{/if}
