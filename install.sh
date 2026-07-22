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
    RC_FILE="$HOME/.zshrc"
fi
touch "$RC_FILE"
echo " [2/4] Detecting shell configuration..."
echo "       Target: ${RC_FILE}"

backup_rc_file() {
    local backup_file="${RC_FILE}.gpumdkit.bak.$(date +%Y%m%d_%H%M%S)"
    cp "$RC_FILE" "$backup_file"
    echo "       Backup created: ${backup_file}"
}

remove_old_gpumdkit_config() {
    local tmp_file
    tmp_file=$(mktemp)

    # Remove the managed GPUMDkit block first.
    awk '
        /^########### GPUMDkit Configuration ###########$/ { in_block=1; next }
        /^##############################################$/ && in_block { in_block=0; next }
        !in_block { print }
    ' "$RC_FILE" > "$tmp_file"

    # Remove older single-line GPUMDkit entries if they were not inside the block.
    grep -v -E '(^export GPUMDkit_path=|^export PATH=\$\{GPUMDkit_path\}:\$\{PATH\}$|^source \$\{GPUMDkit_path\}/Scripts/utils/completion\.sh$)' "$tmp_file" > "${tmp_file}.clean"
    mv "${tmp_file}.clean" "$RC_FILE"
    rm -f "$tmp_file"
}

write_gpumdkit_config() {
    echo "       Adding environment variables to ${RC_FILE}"
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
}

# 3. Write environment variables
if grep -q "export GPUMDkit_path=" "$RC_FILE"; then
    current_paths=$(grep "export GPUMDkit_path=" "$RC_FILE" | sed 's/^[[:space:]]*export GPUMDkit_path=//')
    echo "       Existing GPUMDkit configuration found."
    echo "       Existing path(s):"
    echo "$current_paths" | sed 's/^/         - /'
    echo "       New path:"
    echo "         - ${INSTALL_DIR}"
    echo ""
    read -r -p "       Replace the existing GPUMDkit configuration with the new path? [y/N]: " replace_config

    if [[ "$replace_config" == [yY]* ]]; then
        backup_rc_file
        remove_old_gpumdkit_config
        write_gpumdkit_config
    else
        echo "       Keeping existing GPUMDkit configuration."
    fi
else
    backup_rc_file
    write_gpumdkit_config
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
