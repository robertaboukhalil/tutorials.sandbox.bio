import { SomeFile } from "../utils.js";

let Bed1 = new SomeFile("a.bed", "default", `
chr1	10	20	a1	1	+
chr1	50	70	a2	2	+
chr1	80	90	a3	2	+`);

let Bed2 = new SomeFile("b.bed", "default", `
chr1	30	40	b1	2	+
chr1	55	65	b2	2	+
chr1	85	120	b3	4	+`);


export const Lessons = [
	{
		id: "lesson-loading",
		inputs: [],
	},
	{
		id: "intersect-1",
		title: "bedtools intersect to do xyz",
		description: "that generates overlapping regions between <code>a.bed</code> and <code>b.bed</code>",
		usage: "intersect --help",
		tool: "intersect",
		inputs: [ Bed1, Bed2 ],
		goal: "intersect -a a.bed -b b.bed",
	},
	{
		id: "merge-1",
		title: "bedtools merge to do abc",
		description: "that merges overlapping regions between the two <code>.bed</code> files",
		usage: "merge --help",
		tool: "merge",
		inputs: [ Bed1 ],
		goal: "merge --help",
	}
];
