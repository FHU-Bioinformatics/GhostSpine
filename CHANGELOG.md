# Change Log

## [1.0.0] - 2026-06-29
* Base release

## [1.1.0] - 2026-06-30
* Added changelog
* Added generation of T-free regions from a sequence, as well as visualization and summary stats
* Added read length filtration in aggregate mode
* Improved visualization of sequence Data Frame
* Searching for a read by name is now easier
* A read is now only searched for when the submit button is pressed, instead of any time the input changed

## [1.2.0] - 2026-07-06
* Added chart for Uracil makeup in aggregate mode
* Added T-free boolean mask to read visualization Data Frame
* Added Q-score comparison for T-free and T-bearing regions
* T-free visualization will not be shown if no bases fit the T-free criteria
* Added backend support for analyzing bearing and free regions of any base
* Free regions are now defined as the neighboring 3 bases to the target instead of 6
* Added histogram visualization for U-free and U-bearing regions

## [1.2.1] - Unreleased
* Added Q-score analysis for U-Free and U-Bearing regions in aggregate mode
* Added option to apply coloring to a canonical sequence
* Removed tooltips with no function from read search input forms
* Replaced tkinter file selection with crossfiledialog