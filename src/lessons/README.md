### Lesson Format

```javascript
{
   // Use unique ID (not array position) for saving user progress. This is to
   // future proof against adding or moving lessons to a different order.
   id: "unique-id",
   // Title of the lesson shown in the dropdown at the top of the page
   title: "title",
   // A short description of what the user should accomplish in this lesson
   description: "description",
   // The command to show usage documentation for
   usage: "intersect --help",
   // The name of the program that users should use for this lesson
   tool: "intersect",

   // An array of SomeFile objects that denote input files
   inputs: [ Bed1, Bed2 ],

    // Bedtools command that gives the expected answer (executed at lesson load time)
   goal: "intersect -a a.bed -b b.bed",

    // Do not use. App.svelte will fill this in with the user's answer
   // answer: ""
}
```
