#
# This script is part of GPUMDkit.
# Repository: https://github.com/zhyan0603/GPUMDkit
#
# Description:
#     Update GPUMDkit to latest version
#
# Usage:
#     bash update_gpumdkit.sh [arguments]
#
# Author: Zihan YAN
# Contact: yanzihan@westlake.edu.cn
# Last Modified: 2025-12-28
#

function update_gpumdkit(){
# Check if in a Git repository; if not, try to switch to GPUMDkit_path
if ! git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
    if [ -z "$GPUMDkit_path" ]; then
        echo "Error: Not in a Git repository and GPUMDkit_path environment variable is not set."
        exit 1
    fi
    if [ ! -d "$GPUMDkit_path" ]; then
        echo "Error: Directory $GPUMDkit_path does not exist."
        exit 1
    fi
    cd "$GPUMDkit_path" || {
        echo "Error: Failed to change directory to $GPUMDkit_path."
        exit 1
    }
    if ! git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
        echo "Error: $GPUMDkit_path is not a Git repository."
        exit 1
    fi
fi

# Get current branch name
current_branch=$(git rev-parse --abbrev-ref HEAD)
if [ -z "$current_branch" ]; then
    echo "Error: Unable to determine current branch."
    exit 1
fi

# Get local commit hash for current branch
local_commit=$(git rev-parse HEAD)

# Get remote commit hash for current branch
remote_commit=$(git ls-remote https://github.com/zhyan0603/GPUMDkit.git "$current_branch" | awk '{print $1}')

# Check if remote commit was retrieved successfully
if [ -z "$remote_commit" ]; then
    echo "Error: Unable to fetch remote repository information."
    echo "Check network or branch name ($current_branch)."
    exit 1
fi

# Compare local and remote commits
if [ "$local_commit" = "$remote_commit" ]; then
    echo "No updates available: Local $current_branch branch is up to date."
else
    echo "Updates detected, pulling latest code for branch $current_branch..."
    # Update permissions and pull code
    if [ -f "gpumdkit.sh" ]; then
        chmod -x gpumdkit.sh
    else
        echo "Warning: gpumdkit.sh not found, skipping permission change."
    fi

    # Pull latest code
    if git pull origin "$current_branch"; then
        echo "Code successfully updated."
        # Restore executable permission for gpumdkit.sh
        if [ -f "gpumdkit.sh" ]; then
            chmod +x gpumdkit.sh
            echo "Restored executable permission for gpumdkit.sh."
        fi
    else
        echo "Error: Failed to pull code. Check Git configuration or network connection."
        exit 1
    fi
fi
}