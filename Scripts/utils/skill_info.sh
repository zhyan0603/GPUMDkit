#!/bin/bash

# ============================================================
# GPUMDkit agent skill information helper
# Repository: https://github.com/zhyan0603/GPUMDkit
# Author: Zihan YAN (yanzihan@westlake.edu.cn)
# ============================================================

function skill_info_table(){
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " |                                       GPUMDkit Agent Skill                                            |"
    echo " +-------------------------------------------------------------------------------------------------------+"
    echo " Skill source path:"
    echo "   ${GPUMDkit_path}/skills/gpumdkit-skill"
    echo ""
    echo " Available skills:"
    echo "   gpumdkit-skill                 GPUMDkit, GPUMD, NEP, simulation, and post-processing (English)"
    echo "   gpumdkit-skill-zh              GPUMDkit, GPUMD, NEP, simulation, and post-processing (Chinese)"
    echo ""
    echo " Recommended agent behavior:"
    echo "   Ask about unresolved scientific choices; never guess simulation parameters."
    echo "   Do not launch simulations, training, or scheduler jobs without an explicit request."
    echo "   Ask the user whether to install skills globally or only for the current project."
    echo "   Global install is usually better when GPUMDkit is used across multiple projects."
    echo "   Project install is better when the user wants repository-local configuration."
    echo ""
    echo " Global install targets:"
    echo "   Agent Skills:   ~/.agents/skills"
    echo "   OpenCode:       ~/.config/opencode/skills"
    echo "   Claude compat:  ~/.claude/skills"
    echo ""
    echo " Project install targets:"
    echo "   Agent Skills:   .agents/skills"
    echo "   OpenCode:       .opencode/skills"
    echo "   Claude compat:  .claude/skills"
    echo ""
    echo " Example: install globally using the cross-client convention:"
    echo "   target_dir=\"\${HOME}/.agents/skills\""
    echo "   mkdir -p \"\${target_dir}\""
    echo "   ln -s \"\${GPUMDkit_path}/skills/gpumdkit-skill\" \"\${target_dir}/gpumdkit-skill\""
    echo "   ln -s \"\${GPUMDkit_path}/skills/gpumdkit-skill-zh\" \"\${target_dir}/gpumdkit-skill-zh\""
    echo ""
    echo " Example: install for the current project only:"
    echo "   target_dir=\".agents/skills\""
    echo "   mkdir -p \"\${target_dir}\""
    echo "   ln -s \"\${GPUMDkit_path}/skills/gpumdkit-skill\" \"\${target_dir}/gpumdkit-skill\""
    echo "   ln -s \"\${GPUMDkit_path}/skills/gpumdkit-skill-zh\" \"\${target_dir}/gpumdkit-skill-zh\""
    echo ""
    echo " Replace target_dir with another install target listed above if needed."
    echo " Read ${GPUMDkit_path}/skills/README.md for details."
    echo " +-------------------------------------------------------------------------------------------------------+"
}
