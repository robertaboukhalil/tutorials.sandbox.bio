<script>
// Exports
export let beds = [];  // Format: [{ name: "test.bed", contents: "chr1\t123\t456", color: "blue" }, ...]

// Imports
import { min, max } from "d3-array";
import { scaleLinear } from "d3-scale";

// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

let names, xMin, xMax, xScale;
let datasets = [];
let boundaries = [];
let height = 120;
let width = 1000;

let padding = {left: 60, right: 0, top: 10, bottom: 10};
let boxHeight = (height - padding.top - padding.bottom) / 4; // includes the boxGap.
let boxGap = boxHeight * 0.10;

const colorType = {
	input: "#00ABE7",
	goal: "#9393CB",
	correct: "#3BA99C",
	incorrect: "#EA526F"
};

// -----------------------------------------------------------------------------
// Reactive statements:
// -----------------------------------------------------------------------------

// BED file names and intervals
$: beds = beds.filter(bed => bed.contents != "");
$: names = beds.map(bed => bed.name);
$: colors = beds.map(bed => colorType[bed.type || "input"]);
$: datasets = beds.map(bed => parseBed(bed.contents));

// SVG coordinates
$: xMin = min(datasets, bed => min(bed, d => d.start));
$: xMax = max(datasets, bed => max(bed, d => d.end));
$: xScale = scaleLinear().domain([xMin, xMax]).range([padding.left, width - padding.left - padding.right]);
$: boundaries = [...new Set(datasets.map(bed => bed.map(d => [d.start, d.end])).flat(2))];

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
		<!-- Dotted vertical lines at the boundaries of all intervals -->
		{#each boundaries as edge}
			<line
				x1={xScale(edge)} x2={xScale(edge)}
				y1={padding.top} y2={boxHeight * 4}
				stroke-dasharray="3,3"
				style="stroke:gray;stroke-width:1" />
		{/each}
		{#each datasets as row, i}
			<g transform="translate(0, {padding.top + i * boxHeight})">
				<text x="0" y={boxHeight/2 + boxGap} fill="black">{names[i]}</text>
				{#each row as interval}
				<rect
					x={xScale(interval.start)}
					width={xScale(interval.end) - xScale(interval.start)}
					y="0"
					height={boxHeight - boxGap}
					fill={colors[i]}
				/>
				{/each}
			</g>
		{/each}
  </svg>
</div>
