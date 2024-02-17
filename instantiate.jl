using Pkg

basedir = "docs"

for (root, dirs, files) in walkdir(basedir)
    for file in files
        if endswith(file, "Project.toml")
            proj = joinpath(root, file)
            Pkg.activate(Base.current_project(proj))
            Pkg.instantiate()
            Pkg.precompile()
        end
    end
end
