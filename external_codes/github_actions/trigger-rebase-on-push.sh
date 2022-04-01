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
    PR_NUMBER=$(echo "$pr" | jq -r .[0].pr_number)
    PR_USER_LOGIN=$(echo "$pr" | jq -r .[0].pull_request.user.login)

    # Only process if the pr wants to be autorebased, else save some cycles
    if [[ "$BODY" == *"!-> Feel free to automatically rebase this PR. <-!"* ]]; then
        AUTO_MERGE=$(echo "$pr" | jq -r .[0].auto_merge)
        if [[ "$AUTO_MERGE" != "null" ]]; then
            source external_codes/github_actions/auto-rebase.sh
            # edit pr to cause rebase
            UPDATE_PARAMETERS=$(jq --null-input \
            --arg body "${BODY}\r\n"'AUTOMATED CHANGE: Rebase to new base head of '"${HEAD}" \
            '{"body": $body}')

            RESULT=$(curl -X PATCH -H "${AUTH_HEADER}" -H "${API_HEADER}" \
            -d "${UPDATE_PARAMETERS}" \
            "${URI}/repos/$GITHUB_REPOSITORY/pulls/${PR_NUMBER}")
        fi
    fi
done