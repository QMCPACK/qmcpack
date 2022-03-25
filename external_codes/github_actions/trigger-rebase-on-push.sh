set -e


if [[ -z "$GITHUB_TOKEN" ]]; then
	echo "Set the GITHUB_TOKEN env variable."
	exit 1
fi

URI=https://api.github.com
API_HEADER="Accept: application/vnd.github.v3+json"
AUTH_HEADER="Authorization: token $GITHUB_TOKEN"
PARAMETERS='{"state":"open"}'


pr_list=$(curl -X GET -s -H "${AUTH_HEADER}" -H "${API_HEADER}" \
        -d "${PARAMETERS}" \
		"${URI}/repos/$GITHUB_REPOSITORY/pulls")

pr_list=$(echo "${pr_list}" | jq '[.[] | {pr_number: .number, body: .body, head: .base.sha}]')
for pr in "${pr_list[@]}"; do
    BODY=$(echo "$pr" | jq -r .[0].body)
    HEAD=$(echo "$pr" | jq -r .[0].head)
    PULL_NUMBER=$(echo "$pr" | jq -r .[0].pr_number)
    if [[ "$BODY" == *"!-> Feel free to automatically rebase this PR. <-!"* ]]; then
        # edit pr to cause rebase
        UPDATE_PARAMETERS='{"body":'
        UPDATE_PARAMETERS+="\"${BODY//$'\n'/'\n'}"
        UPDATE_PARAMETERS+="\r\n"
        UPDATE_PARAMETERS+="AUTOMATED CHANGE: Rebase to new base head of ${HEAD}\"}"
        RESULT=$(curl -X PATCH -H "${AUTH_HEADER}" -H "${API_HEADER}" \
        -d "${UPDATE_PARAMETERS}" \
		"${URI}/repos/$GITHUB_REPOSITORY/pulls/${PULL_NUMBER}")
        RESULT=$(echo "$RESULT")
    fi
done