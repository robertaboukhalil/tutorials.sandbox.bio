<script>
import { Lessons } from "../lessons/bowtie.js";
import Tutorial from "../components/Tutorial.svelte";
import VizBed from "../components/VizBed.svelte";
</script>

<Tutorial
    let:files={files}
	toolID="bowtie"
	toolName="Bowtie"
	toolAioli="bowtie/align-s/2.3.5.1"
	toolCLI="bowtie"
    postprocess={d => d.split("\n").filter(d => !(d.startsWith("@") || d.trim() == "")).join("\n")}
	lessons={Lessons}>

    <!-- Visualize .sam output -->
    <div class="row mt-2" slot="viz">
        <div class="col-12">
            Some viz here
        </div>
    </div>
</Tutorial>
