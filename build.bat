@echo off
call echo building project...
call browserify script.js -o bundle.js --debug
call echo browserify complete...
call terser bundle.js --compress -o bundle.min.js
call echo minify complete...
call echo build complete!