<script>
// TODO:
// - Add ability to show a hint to the user
// - Intro lesson that introduces bedtools and what it is
// - Test in Chrome
// - How does it look on mobile?

import { Lessons } from "./lessons/bedtools.js";
import Tutorial from "./components/Tutorial.svelte";
</script>

<Tutorial
	toolID="bedtools"
	toolName="Bedtools"
	toolAioli="bedtools/2.29.2"
	toolCLI="bedtools"
	lessons={Lessons}
/>
