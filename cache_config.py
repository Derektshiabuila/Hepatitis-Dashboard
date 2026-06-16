from flask_caching import Cache

# Define a shared Cache instance to prevent circular imports in dynamic Dash pages
cache = Cache()
