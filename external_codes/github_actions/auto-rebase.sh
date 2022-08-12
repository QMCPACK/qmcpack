#!/bin/bash

# MIT License

# Copyright (c) 2019 cirrus-actions

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

set -e

if [ -z "$PR_NUMBER" ]; then
	PR_NUMBER=$(jq -r ".pull_request.number" "$GITHUB_EVENT_PATH")
	if [[ "$PR_NUMBER" == "null" ]]; then
		PR_NUMBER=$(jq -r ".issue.number" "$GITHUB_EVENT_PATH")
	fi
	if [[ "$PR_NUMBER" == "null" ]]; then
		echo "Failed to determine PR Number."
		exit 1
	fi
fi

echo "Collecting information about PR #$PR_NUMBER of $GITHUB_REPOSITORY..."

if [[ -z "$GITHUB_TOKEN" ]]; then
	echo "Set the GITHUB_TOKEN env variable."
	exit 1
fi

URI=https://api.github.com
API_HEADER="Accept: application/vnd.github.v3+json"
AUTH_HEADER="Authorization: token $GITHUB_TOKEN"

MAX_RETRIES=${MAX_RETRIES:-6}
RETRY_INTERVAL=${RETRY_INTERVAL:-10}
REBASEABLE=""
pr_resp=""
for ((i = 0 ; i < $MAX_RETRIES ; i++)); do
	pr_resp=$(curl -X GET -s -H "${AUTH_HEADER}" -H "${API_HEADER}" \
		"${URI}/repos/$GITHUB_REPOSITORY/pulls/$PR_NUMBER")
	REBASEABLE=$(echo "$pr_resp" | jq -r .rebaseable)
	if [[ "$REBASEABLE" == "null" ]]; then
		echo "The PR is not ready to rebase, retry after $RETRY_INTERVAL seconds"
		sleep $RETRY_INTERVAL
		continue
	else
		break
	fi
done

if [[ "$REBASEABLE" != "true" ]] ; then
	echo "GitHub doesn't think that the PR is rebaseable!"
	exit 1
fi

BASE_REPO=$(echo "$pr_resp" | jq -r .base.repo.full_name)
BASE_BRANCH=$(echo "$pr_resp" | jq -r .base.ref)

USER_LOGIN=$(jq -r ".comment.user.login" "$GITHUB_EVENT_PATH")
          
if [[ "$USER_LOGIN" == "null" ]]; then
	USER_LOGIN=$(jq -r ".pull_request.user.login" "$GITHUB_EVENT_PATH")
fi

# EDIT FROM ORIGNAL SCRIPT: Add case for if this is triggered off a push request
if [[ "$USER_LOGIN" == "null" ]]; then
	USER_LOGIN=$PR_USER_LOGIN
fi

user_resp=$(curl -X GET -s -H "${AUTH_HEADER}" -H "${API_HEADER}" \
	"${URI}/users/${USER_LOGIN}")

USER_NAME=$(echo "$user_resp" | jq -r ".name")
if [[ "$USER_NAME" == "null" ]]; then
	USER_NAME=$USER_LOGIN
fi
USER_NAME="${USER_NAME} (Rebase PR Action)"

USER_EMAIL=$(echo "$user_resp" | jq -r ".email")
if [[ "$USER_EMAIL" == "null" ]]; then
	USER_EMAIL="$USER_LOGIN@users.noreply.github.com"
fi

if [[ -z "$BASE_BRANCH" ]]; then
	echo "Cannot get base branch information for PR #$PR_NUMBER!"
	exit 1
fi

HEAD_REPO=$(echo "$pr_resp" | jq -r .head.repo.full_name)
HEAD_BRANCH=$(echo "$pr_resp" | jq -r .head.ref)

echo "Base branch for PR #$PR_NUMBER is $BASE_BRANCH"

USER_TOKEN=${USER_LOGIN//-/_}_TOKEN
UNTRIMMED_COMMITTER_TOKEN=${!USER_TOKEN:-$GITHUB_TOKEN}
COMMITTER_TOKEN="$(echo -e "${UNTRIMMED_COMMITTER_TOKEN}" | tr -d '[:space:]')"

git remote set-url origin https://x-access-token:$COMMITTER_TOKEN@github.com/$GITHUB_REPOSITORY.git
# CHANGED FROM ORIGNAL: use bot credentials
git config --global user.email "QMCPACKbot@gmail.com"
git config --global user.name "QMCPACK-Bot"

git remote add fork https://x-access-token:$COMMITTER_TOKEN@github.com/$HEAD_REPO.git

set -o xtrace

# make sure branches are up-to-date
git fetch origin $BASE_BRANCH
git fetch fork $HEAD_BRANCH

# do the rebase
git checkout -b fork/$HEAD_BRANCH fork/$HEAD_BRANCH
git rebase origin/$BASE_BRANCH

# push back
git push --force-with-lease fork fork/$HEAD_BRANCH:$HEAD_BRANCH

# CHANGE FROM ORIGINAL: add empty commit signed by bot
echo "$QMCPACK_BOT_GPG_KEY" >> import.key
gpg --import import.key
rm import.key
git config --global user.signingkey "$QMCPACK_BOT_GPG_SIGNING_KEY"  
git commit --allow-empty -S -m "Rebased Signed Off by QMCPACK-Bot"
git push fork fork/$HEAD_BRANCH:$HEAD_BRANCH

COMMIT_SHA=$(git rev-parse --verify HEAD)

# EDIT TO ORIGINAL SCRIPT:  Reapprove PR after push
# NOTE: Submits a review appoving the SHA it just sent, 
# so if someone tried to commit something on top and "steal" this review,
# it shouldn't work because it would have a different SHA
REVIEW_URI="${URI}/repos/$GITHUB_REPOSITORY/pulls/$PR_NUMBER/reviews"
UPDATE_PARAMETERS=$(jq --null-input \
            --arg commit_sha "${COMMIT_SHA}" \
            '{"commit_id" : $commit_sha, "body":"AUTOMATED REVIEW: Reapprove after rebase", "event":"APPROVE"}')


RESULT=$(curl -X POST -H "${AUTH_HEADER}" -H "${API_HEADER}" \
            -d "${UPDATE_PARAMETERS}" \
            "${REVIEW_URI}")