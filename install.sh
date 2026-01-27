#!/bin/bash

# ============================================================
# GPUMDkit Installation Script
# ============================================================

# 1. Get the absolute path of GPUMDkit
# ------------------------------------------------------------
INSTALL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
echo "Detected GPUMDkit directory: ${INSTALL_DIR}"

# 2. Determine the Shell configuration file
# ------------------------------------------------------------
# Defaults to .bashrc; switches to .zshrc if the current shell is zsh
RC_FILE="$HOME/.bashrc"
if [[ "$SHELL" == *"zsh"* ]]; then
    if [ -f "$HOME/.zshrc" ]; then
        RC_FILE="$HOME/.zshrc"
    fi
fi
echo "Target configuration file: ${RC_FILE}"

# 3. Write environment variables
# ------------------------------------------------------------
# Check if the configuration already exists to prevent duplicate entries
if grep -q "export GPUMDkit_path=" "$RC_FILE"; then
    echo -e "\033[33mWarning: GPUMDkit configuration already exists in ${RC_FILE}.\033[0m"
    echo "Skipping environment variable setup to avoid duplicates."
else
    echo "Adding environment variables to ${RC_FILE}..."
    
    # Append the configuration block
    {
        echo ""
        echo "# --- GPUMDkit Configuration ---"
        echo "export GPUMDkit_path=\"${INSTALL_DIR}\""
        echo "export PATH=\"\${GPUMDkit_path}:\${PATH}\""
        
        # Add tab completion support if the script exists
        if [ -f "${INSTALL_DIR}/Scripts/utils/completion.sh" ]; then
            echo "source \"\${GPUMDkit_path}/Scripts/utils/completion.sh\""
        fi
        echo "# ------------------------------"
    } >> "$RC_FILE"
    
    echo -e "\033[32mSuccess: Environment variables added.\033[0m"
fi

# 4. Set executable permissions
# ------------------------------------------------------------
if [ -f "${INSTALL_DIR}/gpumdkit.sh" ]; then
    chmod +x "${INSTALL_DIR}/gpumdkit.sh"
    echo "Added executable permission to gpumdkit.sh"
else
    echo -e "\033[31mError: gpumdkit.sh not found in ${INSTALL_DIR}!\033[0m"
fi

# 5. Final instructions
# ------------------------------------------------------------
echo "------------------------------------------------------------"
echo -e "\033[32mInstallation Complete!\033[0m"
echo ""
echo "To start using GPUMDkit, please run the following command to refresh your shell:"
echo -e "\033[1;33m    source ${RC_FILE}\033[0m"
echo "------------------------------------------------------------"