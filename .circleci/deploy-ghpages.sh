#!/bin/bash
set -eu

SOURCE_DIR="${1:-html}"
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory $SOURCE_DIR does not exist"
    exit 1
fi

git config --global user.email "${GH_EMAIL}"
git config --global user.name "${GH_NAME}"

TEMP_DIR=$(mktemp -d)
trap "rm -rf '$TEMP_DIR'" EXIT
cd "$TEMP_DIR"

git init -b gh-pages

cp -R "${CIRCLE_WORKING_DIRECTORY}/${SOURCE_DIR}"/* .

touch .nojekyll

git add -A
git commit -m "Deploy documentation to GitHub Pages [ci skip]

Built from commit ${CIRCLE_SHA1} on branch ${CIRCLE_BRANCH}"

git remote add origin "https://${GH_TOKEN}@github.com/${CIRCLE_PROJECT_USERNAME}/${CIRCLE_PROJECT_REPONAME}.git"

echo "Pushing to gh-pages branch..."
git push --force --quiet origin gh-pages > /dev/null 2>&1

echo "Documentation deployed successfully to gh-pages"
