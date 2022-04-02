#!/bin/bash

rsync -ruz --progress turso:/proj/wolli/results/ ./results
rsync -ruz --progress turso:/proj/wolli/pbo-instances/ ./instances
