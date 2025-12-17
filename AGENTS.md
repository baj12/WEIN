# AGENTS.md

This file provides guidance to agents when working with code in this repository.


## 🚨 MANDATORY ERROR PROTOCOL 🚨

**READ THIS FIRST: LLMs naturally jump to solutions. You must force yourself to investigate first.**

### When You Encounter An Error

#### ❌ DO NOT:
- Immediately propose a fix
- Add guard clauses (NULL checks, length checks) without investigation
- Assume what values are based on variable names
- Pattern-match error messages to common solutions ("NULL error → add NULL check")

#### ✅ DO (IN THIS ORDER):

**1. State what you DON'T know**
```r
cat("===== ERROR INVESTIGATION =====\n", file = stderr())
cat("ERROR:", error_message, "\n", file = stderr())
cat("UNKNOWN: What is the actual value of problematic_var?\n", file = stderr())
cat("UNKNOWN: Where was it created?\n", file = stderr())
cat("UNKNOWN: What should it be at this point?\n", file = stderr())
```

**2. Trace the call chain BACKWARD**
```r
# Don't assume - TRACE from error point to origin:
# failing_function() line 123
#   ↑ called from: parent_function() line 456
#     ↑ called from: entry_point() line 789
#
# What arguments were passed at each step?
# When did the value become invalid?
```

**3. Insert diagnostic cat() statements**
```r
# At function ENTRY:
cat("=== ENTERING function_name ===\n", file = stderr())
cat("param class:", class(param), "\n", file = stderr())
cat("param is.null:", is.null(param), "\n", file = stderr())
if (!is.null(param)) {
  str(param, max.level = 1, file = stderr())
}
```

**4. Save state BEFORE the failing operation**
```r
# RIGHT BEFORE the line that crashes:
cat("About to execute failing operation\n", file = stderr())
tFile <- tempfile(pattern = "debug_before_crash_", fileext = ".RData")
save.image(tFile)
cat("State saved to:", tFile, "\n", file = stderr())

# Now let it fail - you can load this state and inspect interactively
```

**5. Run with diagnostics and observe**
```r
# Run the code with your cat() statements
# Compare actual output to expectations
# ONLY THEN form hypothesis
```

**6. Propose fix based on evidence**
```r
cat("===== HYPOTHESIS (based on evidence) =====\n", file = stderr())
cat("The value is NULL because [specific reason from diagnostics]\n", file = stderr())
cat("This happens when [specific condition]\n", file = stderr())
cat("Fix: [addresses root cause, not symptom]\n", file = stderr())
```

---

### 📋 Investigation Checklist

**Complete ALL items before proposing any fix:**

- [ ] I have stated what I DON'T know about the error
- [ ] I have traced the call chain backward from error point
- [ ] I have identified WHERE the problematic value originates
- [ ] I have inserted `cat()` diagnostics at key points
- [ ] I have run the code WITH diagnostics and seen actual output
- [ ] I have saved state with `save.image()` before the failing operation
- [ ] I understand the ROOT CAUSE (not just the symptom)
- [ ] My proposed fix addresses root cause, not just adds guards

**If you cannot check all boxes, you are NOT ready to propose a fix.**

---

### ⚠️ Pattern Recognition Trap

LLMs are trained on millions of code examples. When you see an error, instinct triggers pattern-matching:

| Error Message | Instinct Says | **WRONG** |
|--------------|---------------|-----------|
| "NULL" | Add `if (is.null(x))` check | Why is it NULL? |
| "subscript out of bounds" | Add length check | Why is index invalid? |
| "file not found" | Add `file.exists()` | Why is path wrong? |

**These are GUARDS, not FIXES.**

#### When Guards Are Appropriate:
1. You've investigated and found the error is due to **invalid user input**
2. You want to provide a **better error message**
3. The function's **contract explicitly allows** NULL/invalid inputs

#### When Guards Are NOT Appropriate:
1. The error indicates **a bug in YOUR code**
2. The value **SHOULD be valid** but isn't (investigate why!)
3. You **haven't traced** where the invalid value came from

#### Resist The Trap:

When you catch yourself thinking *"I'll just add a check for..."*, **STOP** and ask:

1. **WHY is this value invalid?**
   - Was it created incorrectly?
   - Was it corrupted somewhere?
   - Is the function being called incorrectly?

2. **SHOULD this value be valid here?**
   - What does the function contract say?
   - What do the callers expect?
   - Is this an impossible state that indicates a bug?

3. **WHERE did this value come from?**
   - Trace it backward through the call chain
   - Find where it was created/modified
   - Check if THAT location is the bug

---

### Real Example: NULL Pointer Error

#### ❌ WRONG Approach (Don't Do This)

**Error received:**
```
Error in gs@pointer : no applicable method for '@' applied to an object of class "NULL"
Called from: recompute(gs, pop_path) at line 1687
```

**BAD Response:**
> "Looking at the `compute_boolean_gates` function (lines 1672-1697), the problem is on line 1687 where it's calling `recompute(gs, pop_path)`. The error indicates that `gs` is NULL. Let me add a NULL check at the beginning of the function."

```r
compute_boolean_gates <- function(gs, skip_faulty = FALSE) {
  # WRONG: Just adding guard without understanding WHY
  if (is.null(gs)) {
    stop("gs cannot be NULL")
  }
  # ... rest of function
}
```

**Why this is WRONG:**
- Doesn't investigate WHY `gs` is NULL
- Doesn't check if `gs` SHOULD be NULL
- Doesn't look at the call site to see what's being passed
- Just moves the error earlier - doesn't fix root cause
- No evidence gathered, pure pattern-matching

---

#### ✅ CORRECT Approach (Do This)

**Error received:**
```
Error in gs@pointer : no applicable method for '@' applied to an object of class "NULL"
Called from: recompute(gs, pop_path) at line 1687
```

**GOOD Response:**
> "I see a NULL pointer error. Before proposing a fix, I need to investigate why `gs@pointer` is NULL."

**Step 1: State unknowns**
```r
cat("===== INVESTIGATION START =====\n", file = stderr())
cat("ERROR: gs@pointer is NULL at line 1687 in compute_boolean_gates()\n", file = stderr())
cat("\nUNKNOWNS:\n", file = stderr())
cat("- What is gs at entry to compute_boolean_gates()?\n", file = stderr())
cat("- What does the caller pass as argument?\n", file = stderr())
cat("- Is gs@pointer NULL from the start, or does something corrupt it?\n", file = stderr())
cat("- Does gs@pointer exist but is NULL, or is gs itself NULL?\n", file = stderr())
```

**Step 2: Trace call chain**
```r
cat("\n===== CALL CHAIN TRACE =====\n", file = stderr())
cat("Error location: compute_boolean_gates() line 1687\n", file = stderr())
cat("  recompute(gs, pop_path)\n", file = stderr())
cat("\nFinding caller...\n", file = stderr())
# Search code: compute_boolean_gates is called from fj11_to_gatingset() line 325
cat("Called from: fj11_to_gatingset() line 325\n", file = stderr())
cat("  compute_boolean_gates(gs, skip_faulty = skip_faulty_gate)\n", file = stderr())
cat("\nQUESTION: What is 'gs' at line 325? Need to trace its creation.\n", file = stderr())
```

**Step 3: Insert diagnostics at call site**
```r
# In fj11_to_gatingset(), ADD BEFORE line 325:
cat("\n=== BEFORE compute_boolean_gates() call (line 325) ===\n", file = stderr())
cat("gs class:", class(gs), "\n", file = stderr())
cat("gs is.null:", is.null(gs), "\n", file = stderr())

# Check if gs has @pointer slot
if (!is.null(gs) && isS4(gs)) {
  cat("gs has @pointer slot:", "@pointer" %in% slotNames(gs), "\n", file = stderr())
  if ("@pointer" %in% slotNames(gs)) {
    cat("gs@pointer value:", gs@pointer, "\n", file = stderr())
    cat("gs@pointer is.null:", is.null(gs@pointer), "\n", file = stderr())
  }
} else {
  cat("gs is not an S4 object or is NULL\n", file = stderr())
}

# Also add diagnostics at ENTRY to compute_boolean_gates()
# In compute_boolean_gates(), ADD AT TOP (after line 1675):
cat("\n=== ENTERING compute_boolean_gates() ===\n", file = stderr())
cat("gs class:", class(gs), "\n", file = stderr())
cat("gs is.null:", is.null(gs), "\n", file = stderr())
if (!is.null(gs) && isS4(gs)) {
  str(gs, max.level = 1, file = stderr())
}
```

**Step 4: Save state before failure**
```r
# In compute_boolean_gates(), ADD BEFORE line 1687:
cat("\n=== About to call recompute() at line 1687 ===\n", file = stderr())
cat("pop_path:", pop_path, "\n", file = stderr())
cat("gs status before recompute:\n", file = stderr())
if (!is.null(gs)) {
  cat("  class:", class(gs), "\n", file = stderr())
  cat("  @pointer:", gs@pointer, "\n", file = stderr())
} else {
  cat("  gs IS NULL\n", file = stderr())
}

tFile <- tempfile(pattern = "debug_recompute_", fileext = ".RData")
save.image(tFile)
cat("DEBUG STATE SAVED TO:", tFile, "\n", file = stderr())
cat("Load with: load('", tFile, "')\n", file = stderr())

# Now let it fail - we have the state saved
```

**Step 5: Run and observe diagnostics**
```r
# Run the code with all cat() statements
# Output will show:
#   - What gs actually is at call site (line 325)
#   - What gs is when entering compute_boolean_gates()
#   - What gs@pointer value is
#   - Exact state saved to .RData file

# Example expected output:
# === BEFORE compute_boolean_gates() call (line 325) ===
# gs class: GatingSet
# gs is.null: FALSE
# gs@pointer: <pointer: 0x12345>
#
# === ENTERING compute_boolean_gates() ===
# gs class: GatingSet
# gs is.null: FALSE
#
# === About to call recompute() ===
# gs@pointer: NULL    <--- AH! It's NULL here but wasn't at call site
```

**Step 6: Form hypothesis from EVIDENCE**
```r
cat("\n===== HYPOTHESIS (based on diagnostics output) =====\n", file = stderr())

# If diagnostics show gs@pointer is NULL at call site:
cat("FINDING: gs@pointer is already NULL when passed to compute_boolean_gates()\n", file = stderr())
cat("CONCLUSION: Problem is EARLIER - wherever gs is created\n", file = stderr())
cat("NEXT: Trace back to gs creation in fj11_to_gatingset()\n", file = stderr())

# OR if diagnostics show gs@pointer is valid at entry but NULL at line 1687:
cat("FINDING: gs@pointer is valid at entry but NULL before recompute()\n", file = stderr())
cat("CONCLUSION: Something between entry and line 1687 corrupts it\n", file = stderr())
cat("NEXT: Check operations on gs between lines 1675-1687\n", file = stderr())

# Load the saved state interactively
cat("\nTo debug interactively:\n", file = stderr())
cat("  load('", tFile, "')\n", file = stderr())
cat("  ls()  # see all variables\n", file = stderr())
cat("  str(gs)  # inspect gs structure\n", file = stderr())
```

**Step 7: Fix root cause (example)**
```r
# If investigation revealed that gs@pointer becomes NULL when gs is subset/modified:
# The fix is NOT to add a NULL check, but to preserve the pointer correctly

# BEFORE (broken):
gs_subset <- gs[samples]  # This might invalidate @pointer

# AFTER (fixed):
gs_subset <- gs[samples]
# Re-initialize pointer if needed
if (is.null(gs_subset@pointer)) {
  gs_subset <- reloadGS(gs_subset)  # or whatever the proper fix is
}
```

---
### Most Critical Rules (Top 5)

1. **🛑 Investigate before fixing** - Complete the investigation checklist for every error
2. **🧪 Test-first** - All tests must pass before commits (`testthat::test_check("CyFj11")`)
3. **📊 Real data primary** - Test with actual FlowJo workspaces, mocks are secondary
4. **💾 Save state** - Use `save.image(tempfile())` before failures, debug from saved state
5. **📝 Scripts only** - Never test inline, always create `.R` files in `tests/manual/`


## Stack
- R package using Bioconductor, Shiny, and Shiny Dashboard
- Dependencies include DESeq2, ggplot2, DT, plotly, and many Bioconductor packages
- Tests use testthat framework

## Build/Test Commands
- Run tests: `Rscript -e "devtools::test()"`
- Check package: `Rscript -e "devtools::check()"`
- Install package: `R CMD INSTALL .` or `devtools::install()`
- Run app: `WEIN()` after loading library
- Run tests with data: Source runTest.R after ensuring testData.RData exists


## Git Workflow - Error-Free Commits Only

### Pre-Commit Checklist

```r
# Run this before EVERY commit
source("tests/run_all_checks.R")  # If exists, or:

# Manual checklist:
testthat::test_check("CyFj11")           # All tests pass
devtools::check()                         # R CMD check clean
lintr::lint_package()                     # No linting errors (or acceptable)
```

### Commit Only When

- [ ] All tests pass with zero failures
- [ ] No new warnings introduced
- [ ] Documentation updated (`devtools::document()` run)
- [ ] NAMESPACE auto-generated correctly
- [ ] Changes tested with real data (not just mocks)

### Git Commands for Clean Workflow

```bash
# Before starting work
git status  # Ensure clean working directory

# After changes, verify first
Rscript -e "testthat::test_check('CyFj11')"

# Only if successful:
git add <files>
git commit -m "Descriptive message with test confirmation"

# If tests fail:
git diff  # Review changes
git restore <files>  # Revert if needed
```

---

## Honest Error Reporting - No Fake Success

### Prohibited Patterns

```r
# ❌ NEVER do this:
tryCatch(
  risky_operation(),
  error = function(e) {
    cat("Operation completed successfully\n")  # LIE
    return(NULL)
  }
)

# ❌ NEVER do this:
if (is.null(result)) {
  cat("No errors found\n")  # NULL doesn't mean success
  return(invisible(NULL))
}

# ❌ NEVER do this:
# Just skip the check
# if (validate_gates(ws)) {  # Commented out because failing
#   stop("Validation failed")
# }
cat("Validation passed\n")
```

### Required Patterns

```r
# ✓ DO this:
result <- risky_operation()
if (is.null(result)) {
  stop("risky_operation() returned NULL - operation failed")
}
cat("Operation succeeded, result has", length(result), "elements\n", file = stderr())

# ✓ DO this:
tryCatch(
  {
    result <- risky_operation()
    cat("Success: processed", nrow(result), "rows\n", file = stderr())
    result
  },
  error = function(e) {
    message("FAILED: ", conditionMessage(e))
    stop(e)  # Re-throw, don't suppress
  }
)

# ✓ DO this:
validation_result <- validate_gates(ws)
if (!validation_result$valid) {
  stop("Validation FAILED: ", paste(validation_result$errors, collapse = "; "))
}
cat("Validation passed:", validation_result$gates_checked, "gates checked\n", file = stderr())
```

### Success Messages Must Include Evidence

```r
# ❌ Vague claim:
cat("Gates extracted successfully\n")

# ✓ Concrete evidence:
cat("Extracted", length(gates), "gates from", length(samples), "samples\n", file = stderr())
cat("Gate types:", paste(unique(sapply(gates, function(g) g$type)), collapse = ", "), "\n", file = stderr())
```

---


## Project Structure
- Modular Shiny app with separate UI and server components
- UI modules in R/*_ui.R files
- Server modules in R/*_server.R files
- Core app functions in R/WEIN.R
- Helper functions in R/helpers.R
- Tests in tests/testthat/

## Critical Patterns
- Uses Shiny modules with namespace pattern (e.g., ui_setup_server/ui_setup_ui)
- Heavy use of reactiveValues for state management
- Bookmarking support via onBookmark/onRestore callbacks
- Extensive use of shinydashboard for UI layout
- All server modules follow pattern: function(input, output, session, values, ...)
- UI components use NS() for namespacing
- Global environment ideal_env used for state persistence
- Extensive use of Bioconductor packages for genomic analysis

## Code Style
- Roxygen2 documentation for all exported functions
- Function names use snake_case
- Module functions follow naming pattern *_server/*_ui
- Heavy use of require() inside functions rather than top-level library()
- Extensive import declarations in NAMESPACE file
- Comments use # for single line, longer comments for sections

## Testing Specifics
- Tests use testthat framework
- Test files named test_*.R in tests/testthat/
- Tests typically create small example datasets with DESeq2::makeExampleDESeqDataSet()
- Integration tests in test_integration.R
- Component tests in test_components.R
- Utility function tests in test_utils.R
- Shiny-specific tests in test_shiny.R

## Gotchas
- Many Bioconductor dependencies must be installed separately
- testData.RData required for some test scripts but not committed to repo
- Global variable ideal_env used for state persistence across sessions
- Extensive use of require() inside functions means dependencies checked at runtime
- Package assumes many org.*.eg.db packages may be installed

## Working with Saved States (WEINState files)

### Loading and Using WEINState Files

WEINState files (e.g., `WEINState_20251214_183511.43207.RData`) contain saved application states that can be used for:
- Debugging specific issues
- Testing functionality with real data
- Reproducing user-reported problems

To load and work with a WEINState file:

```r
# Load the saved state
load("WEINState_20251214_183511.43207.RData")

# The state is typically stored in the r_data object
# Access key components:
dds_obj <- r_data$dds_obj        # DESeq2 dataset object
res_obj <- r_data$res_obj        # Results object
cur_species <- r_data$cur_species # Current species selection
cur_type <- r_data$cur_type      # Current gene identifier type

# Example: Inspect the data
head(res_obj)
dim(dds_obj)
cur_species
cur_type
```

### Pattern for Testing with Saved States

When investigating issues or testing functionality:

1. **Load the state first**:
   ```r
   load("WEINState_20251214_183511.43207.RData")
   ```

2. **Extract required data objects**:
   ```r
   # Always check what's available
   ls()           # List objects in workspace
   names(r_data)   # List components of r_data
   ```

3. **Create minimal reproduction scripts**:
   ```r
   # Pattern for testing specific functions
   library(WEIN)
   
   # Load state
   load("WEINState_20251214_183511.43207.RData")
   
   # Extract data
   test_data <- r_data$some_component
   
   # Test function with real data
   result <- some_wein_function(test_data)
   
   # Examine results
   str(result)
   ```

4. **Save intermediate states for debugging**:
   ```r
   # During investigation, save states at key points
   save.image("debug_state_before_issue.RData")
   
   # Later, load and inspect
   load("debug_state_before_issue.RData")
   # Continue investigation...
   ```

### Common Data Objects in WEINState Files

| Object Name | Description | Typical Usage |
|-------------|-------------|----------------|
| `r_data$dds_obj` | DESeq2 dataset object | Differential expression analysis |
| `r_data$res_obj` | DESeq2 results object | Results exploration and visualization |
| `r_data$cur_species` | Current species selection | Species-specific operations |
| `r_data$cur_type` | Gene identifier type | Gene mapping and annotation |
| `r_data$genelistUP` | Upregulated gene list | Functional enrichment analysis |
| `r_data$genelistDOWN` | Downregulated gene list | Functional enrichment analysis |
| `r_data$genelistUPDOWN` | Combined gene list | Comprehensive analysis |

This pattern ensures consistent and reliable testing with real user data while maintaining the ability to reproduce and debug issues effectively.