<script>
import { Lessons } from "../lessons/bedtools.js";
import Tutorial from "../components/Tutorial.svelte";
import VizBed from "../components/VizBed.svelte";
</script>

<Tutorial
    let:files={files}
	toolID="bedtools"
	toolName="Bedtools"
	toolAioli="bedtools/2.29.2"
	toolCLI="bedtools"
	lessons={Lessons}>

    <!-- Visualize .bed files -->
    <div class="row mt-2" slot="viz">
        <div class="col-12">
            <VizBed beds={files} />
        </div>
    </div>
</Tutorial>
