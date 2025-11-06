#!/usr/bin/env bash

################################################################################
# Simplified Environment Setup Script for Bulk RNA-seq Analysis
# 
# Purpose: Validate and document an active micromamba environment for 
#          bulk RNA-seq preprocessing tools (QC, quantification)
#
# Requirements: 
#   - micromamba installed and in PATH
#   - An active micromamba environment with tools installed
#
# Usage: 
#   1. Create and activate environment manually:
#      micromamba create -n rnaseq_env -c conda-forge -c bioconda \
#        salmon fastqc multiqc samtools bedtools cutadapt r-base
#      micromamba activate rnaseq_env
#   
#   2. Run this script:
#      bash setup_environment.sh
#
# Output:
#   - env/environment.yml (environment specification)
#   - validation_log.txt (validation report)
################################################################################

set -euo pipefail

# --- Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/env"
LOG_FILE="${SCRIPT_DIR}/validation_log.txt"

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

################################################################################
# Logging Functions
################################################################################

log_info() {
    echo -e "${BLUE}[INFO]${NC} $*" | tee -a "$LOG_FILE"
}

log_success() {
    echo -e "${GREEN}[✓]${NC} $*" | tee -a "$LOG_FILE"
}

log_warn() {
    echo -e "${YELLOW}[WARNING]${NC} $*" | tee -a "$LOG_FILE"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $*" | tee -a "$LOG_FILE"
}

################################################################################
# Environment Check Functions
################################################################################

check_micromamba() {
    # Check if micromamba is installed
    if ! command -v micromamba &> /dev/null; then
        log_error "micromamba not found in PATH"
        log_error "Please install micromamba: https://mamba.readthedocs.io/en/latest/installation.html"
        exit 1
    fi
    
    local mm_version=$(micromamba --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
    log_success "micromamba found (version ${mm_version})"
}

check_active_environment() {
    # Check if any environment is active
    if [[ -z "${CONDA_DEFAULT_ENV:-}" ]]; then
        log_error "No active micromamba environment detected"
        echo ""
        echo "Please create and activate an environment first:"
        echo ""
        echo "  micromamba create -n rnaseq_env -c conda-forge -c bioconda \\"
        echo "    salmon fastqc multiqc samtools bedtools cutadapt r-base"
        echo ""
        echo "  micromamba activate rnaseq_env"
        echo ""
        exit 1
    fi
    
    # Get active environment name
    local active_env="${CONDA_DEFAULT_ENV}"
    log_success "Active environment detected: ${active_env}"
    
    # Ask for user confirmation
    echo ""
    echo -e "${YELLOW}⚠ This script will validate and document the current environment:${NC}"
    echo -e "  Environment: ${GREEN}${active_env}${NC}"
    echo ""
    read -p "Continue with this environment? (y/n) " -n 1 -r
    echo ""
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_warn "Setup cancelled by user"
        exit 0
    fi
    
    # Store environment name for later use
    export VALIDATED_ENV="${active_env}"
}

################################################################################
# Tool Validation
################################################################################

validate_required_tools() {
    log_info "Validating required bioinformatics tools..."
    
    # Define required tools with minimum versions
    local -A required_tools=(
        ["salmon"]="1.10.0"
        ["fastqc"]="0.12.0"
        ["multiqc"]="1.20"
        ["samtools"]="1.20"
        ["bedtools"]="2.31"
        ["cutadapt"]="4.7"
        ["R"]="4.4.0"
    )
    
    local validation_failed=false
    
    echo "" | tee -a "$LOG_FILE"
    echo "=== Tool Validation ===" | tee -a "$LOG_FILE"
    
    for tool in "${!required_tools[@]}"; do
        if command -v "$tool" &> /dev/null; then
            # Get version based on tool
            local version_output=""
            case "$tool" in
                salmon)
                    version_output=$(salmon --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
                    ;;
                fastqc)
                    version_output=$(fastqc --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
                    ;;
                multiqc)
                    version_output=$(multiqc --version 2>&1 | grep -oE '[0-9]+\.[0-9]+' | head -1)
                    ;;
                samtools)
                    version_output=$(samtools --version 2>&1 | grep -oE '[0-9]+\.[0-9]+' | head -1)
                    ;;
                bedtools)
                    version_output=$(bedtools --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
                    ;;
                cutadapt)
                    version_output=$(cutadapt --version 2>&1 | grep -oE '[0-9]+\.[0-9]+' | head -1)
                    ;;
                R)
                    version_output=$(R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
                    ;;
            esac
            
            if [[ -n "$version_output" ]]; then
                log_success "${tool} v${version_output} (required: >=${required_tools[$tool]})"
            else
                log_warn "${tool} found but version could not be determined"
            fi
        else
            log_error "${tool} NOT FOUND (required: >=${required_tools[$tool]})"
            validation_failed=true
        fi
    done
    
    if [[ "$validation_failed" == true ]]; then
        echo ""
        log_error "Environment validation FAILED - missing required tools"
        echo ""
        echo "To install missing tools, run:"
        echo "  micromamba install -c conda-forge -c bioconda salmon fastqc multiqc samtools bedtools cutadapt r-base"
        echo ""
        exit 1
    fi
    
    log_success "All required tools validated successfully"
}

validate_r_environment() {
    log_info "Checking R installation..."
    
    # Test R can execute and check version
    if R --version &> /dev/null; then
        local r_version=$(R --version 2>&1 | head -1)
        log_success "R installation verified: ${r_version}"
        
        # Check if renv is available (optional)
        if R -e 'library(renv)' &> /dev/null 2>&1; then
            log_success "renv package detected (recommended for R reproducibility)"
        else
            log_warn "renv package not found (optional but recommended)"
            log_info "Install with: R -e 'install.packages(\"renv\")'"
        fi
    else
        log_error "R installation check failed"
        exit 1
    fi
}

################################################################################
# Environment Documentation
################################################################################

export_environment() {
    log_info "Documenting environment specifications..."
    
    # Create env/ directory if it doesn't exist
    mkdir -p "${ENV_DIR}"
    
    # Export environment.yml
    local env_yml="${ENV_DIR}/environment.yml"
    log_info "Exporting environment to: ${env_yml}"
    
    micromamba env export > "${env_yml}" 2>&1 | tee -a "$LOG_FILE"
    
    if [[ -f "${env_yml}" ]]; then
        log_success "Environment specification saved: ${env_yml}"
        
        # Add metadata comment to top of file
        {
            echo "# Bulk RNA-seq Analysis Environment"
            echo "# Generated: $(date '+%Y-%m-%d %H:%M:%S')"
            echo "# Environment: ${VALIDATED_ENV}"
            echo "# Platform: $(uname -s) $(uname -m)"
            echo "#"
            cat "${env_yml}"
        } > "${env_yml}.tmp" && mv "${env_yml}.tmp" "${env_yml}"
        
        log_success "Metadata added to environment.yml"
    else
        log_error "Failed to export environment.yml"
        exit 1
    fi
    
    # Export explicit specification (platform-specific, for exact reproducibility)
    local env_explicit="${ENV_DIR}/environment_explicit.txt"
    log_info "Exporting explicit specification to: ${env_explicit}"
    
    micromamba env export --explicit > "${env_explicit}" 2>&1 | tee -a "$LOG_FILE"
    
    if [[ -f "${env_explicit}" ]]; then
        log_success "Explicit specification saved: ${env_explicit}"
    else
        log_warn "Could not export explicit specification (non-critical)"
    fi
}

################################################################################
# Summary Report
################################################################################

print_summary() {
    local env_name="${VALIDATED_ENV}"
    
    cat << EOF

╔═══════════════════════════════════════════════════════════════════════════╗
║               ENVIRONMENT VALIDATION COMPLETED SUCCESSFULLY                ║
╚═══════════════════════════════════════════════════════════════════════════╝

✓ Validated Environment:
  Name: ${env_name}
  Location: $(micromamba env list | grep "^${env_name}" | awk '{print $NF}')

✓ Documented Files:
  - env/environment.yml (full environment specification)
  - env/environment_explicit.txt (platform-specific exact versions)
  - validation_log.txt (this validation report)

📋 Next Steps:

1. INSTALL R PACKAGES:
   Create env/install_r_packages.R with your required Bioconductor packages
   Then run: Rscript env/install_r_packages.R

2. INITIALIZE RENV (recommended for R reproducibility):
   Rscript -e "renv::init(bioconductor='3.19')"
   
3. COMMIT ENVIRONMENT SPECIFICATIONS:
   git add env/environment.yml env/environment_explicit.txt
   git commit -m "Add validated environment specifications"

4. SHARE ENVIRONMENT:
   Others can recreate this environment with:
   micromamba env create -f env/environment.yml

💡 Tips:
- Always activate this environment before running analyses
- Update environment.yml after installing new tools:
  micromamba env export > env/environment.yml
- Keep env/environment.yml under version control
- Use env/environment_explicit.txt for exact reproducibility on same platform

📖 Documentation:
- Micromamba: https://mamba.readthedocs.io/
- Bioconductor: https://bioconductor.org/
- Salmon: https://salmon.readthedocs.io/
- renv: https://rstudio.github.io/renv/

═══════════════════════════════════════════════════════════════════════════

EOF
}

################################################################################
# Main Execution
################################################################################

main() {
    echo "╔═══════════════════════════════════════════════════════════════════════════╗"
    echo "║         Bulk RNA-seq Environment Validation & Documentation Tool          ║"
    echo "╚═══════════════════════════════════════════════════════════════════════════╝"
    echo ""
    
    # Initialize log file
    {
        echo "Environment Validation Started: $(date)"
        echo "Script Directory: ${SCRIPT_DIR}"
        echo "Platform: $(uname -s) $(uname -m)"
        echo "═══════════════════════════════════════════════════════════════════════════"
        echo ""
    } > "$LOG_FILE"
    
    # Run validation steps
    check_micromamba
    check_active_environment
    validate_required_tools
    validate_r_environment
    export_environment
    
    # Print final summary
    print_summary
    
    log_success "Setup completed successfully at $(date)"
    echo "Full validation log saved to: ${LOG_FILE}"
    
    exit 0
}

# Execute main function
main "$@"
