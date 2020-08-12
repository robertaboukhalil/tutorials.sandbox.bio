import { SomeFile } from "../utils.js";

// -----------------------------------------------------------------------------
// BED Files
// -----------------------------------------------------------------------------

let Bed1 = new SomeFile("a.bed", "default", `
chr1	10	20	a1	1	+
chr1	50	70	a2	2	+
chr1	80	90	a3	2	+`);

let Bed2 = new SomeFile("b.bed", "default", `
chr1	30	40	b1	2	+
chr1	55	65	b2	2	+
chr1	85	120	b3	4	+`);


// -----------------------------------------------------------------------------
// Lessons
// -----------------------------------------------------------------------------
// Format:
// {
//    // Use unique ID (not array position) for saving user progress. This is to
//    // future proof against adding or moving lessons to a different order.
//    id: "unique-id",
//    // Title of the lesson shown in the dropdown at the top of the page
//    title: "title",
//    // A short description of what the user should accomplish in this lesson
//    description: "description",
//    // The bedtools command to show usage documentation for
//    usage: "intersect --help",
//    // The name of the bedtools program that users should use for this lesson
//    tool: "intersect",
//    // An array of SomeFile objects that denote input files
//    inputs: [ Bed1, Bed2 ],
//    // Bedtools command that gives the expected answer (executed at lesson load time)
//    goal: "intersect -a a.bed -b b.bed",
//    // Do not use. App.svelte will fill this in with the user's answer
//    // answer: ""
// }
// -----------------------------------------------------------------------------

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
