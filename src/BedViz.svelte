<script>
// Exports
export let beds = [];  // Format: [{ name: "test.bed", contents: "chr1\t123\t456" }, ...]

// Imports
import { min, max } from "d3-array";
import { scaleLinear } from "d3-scale";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

let names, datasets, xMin, xMax, xScale;
let height = 400;
let width = 800;


// -----------------------------------------------------------------------------
// Reactive statements:
// -----------------------------------------------------------------------------

// BED file names and intervals
$: names = beds.map(bed => bed.name);
$: datasets = beds.map(bed => parseBed(bed.contents));

// SVG coordinates
$: xMin = min(datasets, bed => min(bed, d => d.start));
$: xMax = max(datasets, bed => max(bed, d => d.start))
$: xScale = scaleLinear().domain([xMin, xMax]).range([0, width]);


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Parse a .bed file
function parseBed(bed) {
	if(bed == null)
		return [];
	let lines = bed.split("\n");
	return lines.filter(l => l != "").map(parseLine);
}

// Extract info from a .bed file line
function parseLine(line) {
	let l = line.split("\t");
	return {
		chr: l[0],
		start: Number(l[1]),
		end: Number(l[2]),
		id: l[3]
	};
}
</script>

<div height={height} width={width}>
  <svg height={height} width={width}>
	{#each datasets as row, i}
		<g transform="translate(70, {i * 30})">
			<text x="-70" y="20" fill="black">{names[i]}</text>
			{#each row as interval}
			<rect
				x={xScale(interval.start)}
				width={xScale(interval.end) - xScale(interval.start)}
				y="0"
				height="25"
				fill="blue"
			/>
			{/each}
		</g>
	{/each}
  </svg>
</div>
