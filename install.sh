#!/bin/bash

echo "======================================================"
echo "  GPUMDkit Installation"
echo "======================================================"

# 1. Get the absolute path of GPUMDkit
INSTALL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
echo " [1/4] Detecting GPUMDkit directory..."
echo "       ${INSTALL_DIR}"

# 2. Determine the Shell configuration file
RC_FILE="$HOME/.bashrc"
if [[ "$SHELL" == *"zsh"* ]]; then
    if [ -f "$HOME/.zshrc" ]; then
        RC_FILE="$HOME/.zshrc"
    fi
fi
echo " [2/4] Detecting shell configuration..."
echo "       Target: ${RC_FILE}"

# 3. Write environment variables
# Check if the configuration already exists to prevent duplicate entries
if grep -q "export GPUMDkit_path=" "$RC_FILE"; then
    echo "       Warning: GPUMDkit configuration already exists."
    echo "       Skipping environment variable setup to avoid duplicates."
else
    echo "       Adding environment variables to ${RC_FILE}"

    # Append the configuration block
    {
        echo ""
        echo "########### GPUMDkit Configuration ###########"
        echo "export GPUMDkit_path=${INSTALL_DIR}"
        echo "export PATH=\${GPUMDkit_path}:\${PATH}"

        # Add tab completion support if the script exists
        if [ -f "${INSTALL_DIR}/Scripts/utils/completion.sh" ]; then
            echo "source \${GPUMDkit_path}/Scripts/utils/completion.sh"
        fi
        echo "##############################################"
    } >> "$RC_FILE"

    echo "       Success: Environment variables added."
fi

# 4. Set executable permissions
echo " [3/4] Setting executable permissions..."
if [ -f "${INSTALL_DIR}/gpumdkit.sh" ]; then
    chmod +x "${INSTALL_DIR}/gpumdkit.sh"
    echo "       Added executable permission to gpumdkit.sh"
else
    echo "       Error: gpumdkit.sh not found in ${INSTALL_DIR}!"
fi

# 5. Make variables available in current shell
echo " [4/4] Loading environment..."
source "${RC_FILE}"
echo ""
echo "======================================================"
echo "  Installation Complete!  GPUMDkit is ready to use."
echo "======================================================"
echo ""
echo "  Usage:"
echo "    gpumdkit.sh        Interactive mode"
echo "    gpumdkit.sh -h     Show help"
echo "    gpumdkit.sh -<opt> Command-line mode"
echo ""
echo "======================================================"
