cacert.pem is the Mozilla list of trusted certificate authorities, taken from the Python package "certifi"
(version 2026.10.23, downloaded 2026-11-02). MR Studio uses it only as a last resort, so that https access to OpenGWAS
works on a computer where the bundled R cannot find any certificate file of its own. Refresh it from time to time.
