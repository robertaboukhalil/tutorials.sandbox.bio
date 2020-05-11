<script>
import { scaleLinear } from 'd3-scale';
import { min, max } from 'd3-array';

// Received properties:
export let beds = "";
export let names = "";

// Internal globals:
let Datasets = [];
let XScale;
let Width = 800;

// Reactive statements:
$: show(beds);

function parse_line(line) {
  let l = line.split('\t');
  return {chr: l[0], start: Number(l[1]), end: Number(l[2]), id: l[3]}
}
function parse_bed(bed) {
  let lines = bed.split('\n');
  return lines.filter(l => l != "").map(parse_line);
}
function show(beds) {
  Datasets = beds.filter(bed => bed != "").map(parse_bed);
  console.log("Datasets:", Datasets);
  let xmin = min(Datasets, bed => min(bed, d => d.start))
  let xmax = max(Datasets, bed => max(bed, d => d.start))
  XScale = scaleLinear()
    .domain([xmin, xmax])
    .range([0, Width]);
}
</script>

<div height="400" width="{Width}">
  <svg height="400" width="{Width}">
    {#each Datasets as row, i}
      <g transform="translate(70, {i* 30})">
        <text x="-70" y="20" fill="black">{names[i]}</text>
        {#each row as interval}
          <rect
            x="{XScale(interval.start)}"
            width="{XScale(interval.end) - XScale(interval.start)}"
            y="0"
            height="25"
            fill="blue" />
        {/each}
      </g>
    {/each}
  </svg>
</div>
