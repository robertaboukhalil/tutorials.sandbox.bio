# tutorials.sandbox.bio

## Development

> Make sure you have Node.js version >6 (we currently use v13), otherwise, you'll get the error `Unexpected token {`. Use `node --version` to check which version you have. 
> 
> If you need to upgrade Node.js: `npm cache clean -f && npm install -g n && n stable` (may need `sudo`)

Install dependencies and launch dev mode (automatically compiles + refreshes page on change):

```bash
npm install
npm run dev
```

Then open http://localhost:5000 in your browser.

## Modify Lessons

To update the lessons, modify the JSON in [lessons.js](https://github.com/robertaboukhalil/tutorials.sandbox.bio/tree/master/src/lessons).

## Deploy

For the final build (before deployment), use `npm run build`.
