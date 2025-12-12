# Setup file for shinytest2
#
# This file is executed before any test files are run, and can be used to
# configure shinytest2 options or load helper functions.

library(shinytest2)

# Set default timeout for tests
options(shinytest2.default_timeout = 30000)

# Set default screenshot directory
options(shinytest2.screenshot_dir = "test-figs")

# Other shinytest2 configuration can go here