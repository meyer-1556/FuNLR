# Publish the portable alpha through GitHub in a browser

These steps publish FuNLR **0.5.0a2** at [github.com/meyer-1556/FuNLR](https://github.com/meyer-1556/FuNLR), with the alpha README on the repository landing page and a versioned prerelease. They require repository-owner permissions for visibility changes, but no command-line Git.

Use the companion web-upload archive for source files. Its numbered batches contain no more than 90 files each. A ZIP uploaded into the repository remains a ZIP: GitHub does not extract it into source files. Wheel and source-distribution files belong on the release page, not inside the source tree.

## 1. Review what will become public

The prepared source snapshot is separate from the existing repository's history. Before making a private repository public, review its other branches, tags, older commits, issues, pull requests, release attachments, and Actions logs/artifacts for material that should remain private. Uploading a clean version or deleting a file in a new commit does not remove earlier copies from Git history. A public repository cannot contain a private branch.

Changing visibility exposes the repository's code and Actions history, not just the selected alpha branch. Follow [GitHub's visibility guidance](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/managing-repository-settings/setting-repository-visibility). Keep internal assessment documents and working datasets outside the upload batches.

If the existing history must remain private, retain that repository as private and create a separate public repository containing only the reviewed source snapshot. That gives the project a different URL. Preserving the exact original URL would require a deliberate owner-managed repository rename/replacement plan; do not rename or delete the existing repository merely to follow this guide. Check links, integrations and repository contents before choosing that alternative.

## 2. Create a review branch

1. Sign in and open [FuNLR](https://github.com/meyer-1556/FuNLR).
2. Select the repository's actual `main` branch. This establishes the real repository ancestry for the update.
3. In the branch selector, create `alpha/portable-python` from `main`. If that name already exists, use a new name such as `alpha/portable-python-0.5.0a2`.
4. Confirm the new branch is selected before uploading. Keep all upload commits on that branch until the complete source tree has been reviewed.

## 3. Upload every source batch

1. Extract the web-upload ZIP locally. Locate its numbered batch folders and upload manifest.
2. In Finder, press **Command–Shift–period** to reveal hidden files. Include `.github/`, `.gitignore` and `.dockerignore` wherever they occur in the batches.
3. At the repository root on the review branch, choose **Add file → Upload files**. Drag in the **contents** of the first batch folder.
4. Check that preview paths begin with `src/`, `docs/`, or the intended root filename. They must not contain an enclosing release or batch-folder name.
5. Commit to the review branch with a message such as `Add portable alpha source (batch 1)`. Return to the repository root and repeat for every batch. Shared folder names across batches are expected; the file paths differ.
6. Compare the resulting file list against the supplied upload manifest. Uploads replace matching paths but **do not delete obsolete files**. Review old files absent from the snapshot, delete those being retired through their file/directory menus, and commit those deletions. Preserve unrelated files intentionally retained by the project. These deletions affect the current tree, not old commits.

The prepared batches stay below GitHub's limits of 100 files per upload and 25 MiB per file. Browser uploads do not apply `.gitattributes` transformations. See [GitHub's upload documentation](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository).

## 4. Review and merge the alpha to the landing page

1. Open `README.md` and confirm version **0.5.0a2**. Check that `pyproject.toml`, `src/funlr/`, `examples/`, `tests/`, `docs/`, `install.py` and `.github/workflows/` are present at the repository root. Confirm no upload-batch folders or internal working documents appear there.
2. Open **Actions** and inspect the checks for the final, complete upload commit. Intermediate batch commits may fail because their source tree is incomplete. A local test result does not establish that a hosted check passed. Resolve failures on the complete tree before merging; keep the stated validation limits accurate.
3. Open **Pull requests → New pull request**. Choose base `main` and compare the review branch. Use a title such as `Release the portable FuNLR alpha` and review all changed files and deletions.
4. Once the source review and applicable checks are complete, merge the pull request using the available merge control. If branch rules require another review, follow those rules.
5. Return to the repository's default branch and confirm the alpha README is displayed. Keep `main` as the default branch for this route. Check the merged commit's Actions results as well; a branch-only upload would leave the old default-branch README at the main URL.

## 5. Make the existing repository publicly readable

If the repository is private, complete the whole-repository review in Step 1 first. Then open **Settings → General → Danger Zone → Change repository visibility**, choose **Public**, and read and complete GitHub's owner confirmations. If it is already public, no visibility change is needed. This publishes the existing repository at the same URL; it does not sanitize its history. [GitHub documents the visibility controls and their effects](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/managing-repository-settings/setting-repository-visibility).

Open [the repository URL](https://github.com/meyer-1556/FuNLR) in a signed-out or private browser window. Confirm that the source tree, alpha README, license and installation instructions are readable without authentication.

## 6. Create the versioned alpha prerelease

1. Open **Releases → Draft a new release**. In **Choose a tag**, create `v0.5.0a2`; select **Target: main** after verifying it contains the completed, reviewed alpha. If a tag of that name already exists, inspect its commit before using it. Do not attach this release to an older or partial-upload commit.
2. Set the title to **FuNLR 0.5.0a2 — portable alpha**. Paste the contents of [RELEASE_NOTES.md](RELEASE_NOTES.md) into the description. Preview links and retain the distinction between software demonstrations and independent biological validation.
3. Attach the supplied repository ZIP, wheel, source-distribution tarball and matching checksum file. Do not attach internal audits, working datasets or the browser-upload batches. GitHub also generates source ZIP/tar archives from the selected tag.
4. Select **This is a pre-release**. Save a draft while checking the tag target, notes and attachments. When everything is ready, choose **Publish release**. If immutable releases are enabled, attach all intended assets before publishing.

These controls are described in [GitHub's release documentation](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository). Finish by opening the repository and release page while signed out, downloading one attached archive, and confirming that its filename/version match **0.5.0a2**. The public alpha is then discoverable at the repository URL and downloadable as a named prerelease.
