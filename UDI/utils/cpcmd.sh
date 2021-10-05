find mnt -name "results.mat" | while IFS= read -r NAME; do cp -v "$NAME" "./${NAME//\//_}"; done
