# Reference Data

Reference data used by these analyses totals approximately 19 GB, and thus
cannot be hosted on GitHub. Instead, we have hosted them in a Google storage
bucket.  To obtain these data, please [install gsutil](https://cloud.google.com/storage/docs/gsutil_install) and run the following command in this directory (`ref`):

```
gsutil -m cp -r gs://getzlab-passengerhotspots/* .
```
