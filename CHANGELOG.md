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

## [1.1.1] - Unreleased
* Added chart for Uracil makeup in aggregate mode
* Added T-free boolean mask to read visualization Data Frame
* Added Q-score comparison for T-free and T-bearing regions
* T-free visualization will not be shown if no bases fit the T-free criteria
* Added backend support for analysing bearing and free regions of any base
* Free regions are mow defined as the neighboring 3 bases to the target instead of 6
* Added backend histogram visualization for free and bearing regions