#!/bin/bash

echo "======================================================"
echo "  GPUMDkit Installation"
echo "======================================================"

# 1. Get the absolute path of GPUMDkit ('zsh source' has no BASH_SOURCE, falls back to $0)
INSTALL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )
echo " [1/4] Detecting GPUMDkit directory..."
echo "       ${INSTALL_DIR}"

# 2. Determine the Shell configuration file
RC_FILE="$HOME/.bashrc"
if [[ "$SHELL" == *"zsh"* ]]; then
    RC_FILE="$HOME/.zshrc"
elif [[ "$SHELL" != *bash* ]]; then
    echo " [2/4] Detecting shell configuration..."
    echo "       Warning: unsupported shell '${SHELL}' detected."
    echo "       GPUMDkit scripts run under bash; the configuration will be"
    echo "       written to ${RC_FILE}. If you use another shell (e.g. fish),"
    echo "       add the settings to its config manually."
fi
if ! touch "$RC_FILE" 2>/dev/null; then
    echo " Error: cannot write to ${RC_FILE}. Please check the file permissions."
    exit 1
fi
echo " [2/4] Detecting shell configuration..."
echo "       Target: ${RC_FILE}"

backup_rc_file() {
    local backup_file="${RC_FILE}.gpumdkit.bak.$(date +%Y%m%d_%H%M%S)"
    if ! cp "$RC_FILE" "$backup_file"; then
        echo " Error: failed to back up $RC_FILE."
        return 1
    fi
    echo "       Backup created: ${backup_file}"
}

remove_old_gpumdkit_config() {
    local tmp_file grep_status
    if ! tmp_file=$(mktemp); then
        echo " Error: failed to create a temporary file for $RC_FILE."
        return 1
    fi

    # Remove the managed GPUMDkit block first.
    awk '
        /^########### GPUMDkit Configuration ###########$/ { in_block=1; next }
        /^##############################################$/ && in_block { in_block=0; next }
        !in_block { print }
    ' "$RC_FILE" > "$tmp_file" || {
        echo " Error: failed to read $RC_FILE."
        rm -f "$tmp_file"
        return 1
    }

    # Remove older single-line GPUMDkit entries if they were not inside the block.
    grep -v -E '(^export GPUMDkit_path=|^export PATH="?[$]\{GPUMDkit_path\}:[$]\{PATH\}"?$|^source [$]\{GPUMDkit_path\}/Scripts/utils/completion\.sh$)' "$tmp_file" > "${tmp_file}.clean"
    grep_status=$?
    if [ "$grep_status" -gt 1 ]; then
        echo " Error: failed to clean $RC_FILE."
        rm -f "$tmp_file" "${tmp_file}.clean"
        return 1
    fi
    if ! mv "$tmp_file.clean" "$RC_FILE"; then
        echo " Error: failed to update $RC_FILE."
        rm -f "$tmp_file" "$tmp_file.clean"
        return 1
    fi
    rm -f "$tmp_file"
}

write_gpumdkit_config() {
    echo "       Adding environment variables to ${RC_FILE}"
    {
        echo ""
        echo "########### GPUMDkit Configuration ###########"
        echo "export GPUMDkit_path=\"${INSTALL_DIR}\""
        echo "export PATH=\"\${GPUMDkit_path}:\${PATH}\""

        # Add tab completion support if the script exists
        if [ -f "${INSTALL_DIR}/Scripts/utils/completion.sh" ]; then
            echo "source \${GPUMDkit_path}/Scripts/utils/completion.sh"
        fi
        echo "##############################################"
    } >> "$RC_FILE"

    if [ $? -ne 0 ]; then
        echo " Error: failed to write the configuration to ${RC_FILE}."
        exit 1
    fi
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
    # Portable prompt: 'read -p' is bash-only (zsh reads from a coprocess)
    printf '       Replace the existing GPUMDkit configuration with the new path? [y/N]: '
    if ! IFS= read -r replace_config; then
        echo "       No response received; keeping existing GPUMDkit configuration."
        replace_config=""
    fi

    if [[ "$replace_config" == [yY]* ]]; then
        backup_rc_file || exit 1
        remove_old_gpumdkit_config || exit 1
        write_gpumdkit_config
    else
        echo "       Keeping existing GPUMDkit configuration."
    fi
else
    backup_rc_file || exit 1
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
