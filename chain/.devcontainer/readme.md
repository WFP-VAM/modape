# Storage binding to the Dev Container, either:
* Create a symlink in the chain folder (parent of this .devcontainer folder) named "storage", pointing to a folder you wish to mount as /var/storage
OR
* Create a subfolder named "storage" in the chain folder (parent of this .devcontainer folder; use the setup_storage.sh at your convenience)

# Production deployment
The production image is also build from the Dockerfile. Use the compose .yml to run the processing chain