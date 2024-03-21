#!/usr/bin/env bash

set -o errexit
set -o nounset

command -pV nvidia-smi &> /dev/null || {
	>&2 echo "error: unable to find nvidia-smi command"
	exit 1
}

cuda_compiler_version=$(nvidia-smi -q | grep -i '^CUDA Version' | cut -d: -f2 | xargs)

# function to facilitate version comparison; cf. https://stackoverflow.com/a/37939589
version2int () { echo "$@" | awk -F. '{ printf("%d%02d\n", $1, $2); }'; }

declare -a CUDA_CONFIG_ARGS
if [ "${cuda_compiler_version}" != "None" ]; then
    cuda_compiler_version_int=$(version2int "$cuda_compiler_version") 

    ARCHES=()
    if   [ $cuda_compiler_version_int -eq $(version2int "11.0") ]; then
        ARCHES=(35    50    53 60 61 62 70 72 75)
    elif [ $cuda_compiler_version_int -eq $(version2int "11.1") ]; then
        ARCHES=(35    50    53 60 61 62 70 72 75)
    elif [ $cuda_compiler_version_int -eq $(version2int "11.2") ]; then
        ARCHES=(35    50    53 60 61 62 70 72 75 80 86)
    elif [ $cuda_compiler_version_int -eq $(version2int "11.3") ]; then
        ARCHES=(35    50    53 60 61 62 70 72 75 80 86)
    elif [ $cuda_compiler_version_int -eq $(version2int "11.4") ]; then
        ARCHES=(35    50    53 60 61 62 70 72 75 80 86)
    elif [ $cuda_compiler_version_int -eq $(version2int "11.5") ]; then
        ARCHES=(35 37 50 52 53 60 61 62 70 72 75 80 86 87)
    elif [ $cuda_compiler_version_int -eq $(version2int "11.6") ]; then
        ARCHES=(35 37 50 52 53 60 61 62 70 72 75 80 86 87)
    elif [ $cuda_compiler_version_int -eq $(version2int "11.7") ]; then
        ARCHES=(35 37 50 52 53 60 61 62 70 72 75 80 86 87)
    elif [ $cuda_compiler_version_int -eq $(version2int "11.8") ]; then
        ARCHES=(35 37 50 52 53 60 61 62 70 72 75 80 86 87 89 90)
    elif [ $cuda_compiler_version_int -eq $(version2int "12.0") ]; then
        ARCHES=(      50 52 53 60 61 62 70 72 75 80 86 87 89 90 90a)
    elif [ $cuda_compiler_version_int -eq $(version2int "12.1") ]; then
        ARCHES=(      50 52 53 60 61 62 70 72 75 80 86 87 89 90 90a)
    elif [ $cuda_compiler_version_int -eq $(version2int "12.2") ]; then
        ARCHES=(      50 52 53 60 61 62 70 72 75 80 86 87 89 90 90a)
    elif [ $cuda_compiler_version_int -eq $(version2int "12.3") ]; then
        ARCHES=(      50 52 53 60 61 62 70 72 75 80 86 87 89 90 90a)
    fi

    LATEST_ARCH="${ARCHES[-1]}"
    unset "ARCHES[${#ARCHES[@]}-1]"

    for arch in "${ARCHES[@]}"; do
        CMAKE_CUDA_ARCHS="${CMAKE_CUDA_ARCHS+${CMAKE_CUDA_ARCHS};}${arch}-real"
    done

    CMAKE_CUDA_ARCHS="${CMAKE_CUDA_ARCHS+${CMAKE_CUDA_ARCHS};}${LATEST_ARCH}"

    CUDA_CONFIG_ARGS+=(
        "${CMAKE_CUDA_ARCHS}"
    )
fi

echo -n ${CUDA_CONFIG_ARGS+"${CUDA_CONFIG_ARGS[@]}"}
