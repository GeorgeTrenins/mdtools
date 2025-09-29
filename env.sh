#!/bin/bash
# sets environment variables for GLEFit

if [ -n "$BASH_SOURCE" ]; then
    # Bash or zsh (if BASH_SOURCE is defined)
    script="$BASH_SOURCE[0]"
elif [ -n "$ZSH_VERSION" ]; then
    # Zsh
    script="${(%):-%x}"
else
    # Dash or POSIX sh
    script="$0"
fi

ENV_BASE_DIR=$(cd "$(dirname "$(readlink -f "$script")")" && pwd)

echo "Setting up MD tools paths to base folder $ENV_BASE_DIR"

export PATH=$ENV_BASE_DIR/bin:$PATH
export PYTHONPATH="${ENV_BASE_DIR}/src":$PYTHONPATH

unset ENV_BASE_DIR
