# assumes you have set ZENODO_TOKEN (see showyourwork docs)
# DEPOSITION is the final number of the url of a "new version" on zenodo website
# filename hardcoded, file path to be provided as argument

set -xe

DEPOSITION=$1
FILEPATH=$2
FILENAME="MESA_output.tar.gz"

BUCKET=$(curl -H "Accept: application/json" -H "Authorization: Bearer $ZENODO_TOKEN" "https://www.zenodo.org/api/deposit/depositions/$DEPOSITION" | jq --raw-output .links.bucket)

curl --progress-bar -o /dev/null --upload-file $FILEPATH/$FILENAME $BUCKET/$FILENAME?access_token=$ZENODO_TOKEN
