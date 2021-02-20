# Storage binding to the Dev Container, either:
* Create a symlink in the chain folder (parent of this .devcontainer folder) named "storage", pointing to a folder you wish to mount as /var/storage
OR
* Do nothing: a plain subfolder named "storage" will be created in the parent chain folder

# Production deployment
The production image is also build from the Dockerfile. Use the compose .yml to run the processing chain