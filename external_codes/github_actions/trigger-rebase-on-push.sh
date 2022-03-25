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


pr_list | jq '.[] | {pr_number: .number, body: .body, head: .base.sha}' | while read object; do
    BODY=$(echo "$object" | jq -r .body)
    HEAD=$(echo "$object" | jq -r .head)
    PULL_NUMBER=$(echo "$object" | jq -r .pr_number)
    if [[BODY == *"!-> Feel free to automatically rebase this PR. <-!"*]]; then
        # edit pr to cause rebase
        UPDATE_PARAMETERS="{'body':${BODY}" + "\n" + "AUTOMATED CHANGE: Rebase to new base head of ${HEAD}}"
        RESULT = $(curl -X PATCH -H "${AUTH_HEADER}" -H "${API_HEADER}" \
        -d "${PARAMETERS}" \
		"${URI}/repos/$GITHUB_REPOSITORY/pulls/${PULL_NUMBER}")
    fi
done