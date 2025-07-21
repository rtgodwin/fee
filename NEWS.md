# fee 1.0.0 (2025-07-20)

## Initial CRAN Release

This is the first release of the **fee** package. It estimates the first-exposure effect (FEE) in social programs, using count data models estimated using the **oneinfl** package.

### **Features**
- Provides `fee()` for estimating the FEE.
- Provides `dfee()` for estimating the marginal FEE.
- Provides `feeplot()` for plotting the factual and counterfactual distributions.
- Provides `feepred()` for estimating predicted counts under factual and counterfactual distributions.
- Provides documentation.

### **Testing & Compatibility**
- Successfully tested on:
  - ✅ **Windows** (`windows-latest`) via RHUB and `devtools::check_win_devel()`
  - ✅ **Linux** (`ubuntu-latest`) via RHUB
  - ✅ **macOS** (`macos-13`, `macos-15`, and `macos-arm64`) via RHUB
- Passed:
  - ✅ `devtools::check()`
  - ✅ `devtools::check_win_devel()`
  - ✅ All RHUB checks
  - ✅ GitHub Actions CI
- Only NOTE is from `check_win_devel()` indicating a new submission — expected for first-time CRAN packages.
- No known issues.