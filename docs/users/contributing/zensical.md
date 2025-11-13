---
icon: lucide/rocket
tags:
  - contributing
  - documentation
  - user manual
  - zensical  
---

# User's documentation

For generating this user's documentation, we use [zensical](https://zensical.org/docs/).
Zensical is a modern static site generator designed to simplify building and maintaining project documentation.

## Installation

Zensical is written in Rust and Python, and is published as a [Python package].
We  recommend to use a Python _virtual environment_ when installing [with
`pip`][with-pip] or [with `uv`][with-uv]. Both options automatically install all
necessary dependencies alongside Zensical.

[Python package]: https://pypi.org/project/zensical

!!! note "Prerequisites"
    You need to have Python and a Python package manager installed on your
    system before you install Zensical. We recommend you follow the [Python
    Setup and Usage] instructions for your operating system provided on the
    [Python website].  Modern Python distributions include the `pip` package
    manager, so unless you are developing Python software and use `uv`, this is
    the simplest option to install Zensical on your system.

  [with-pip]: #install-with-pip
  [with-uv]: #install-with-uv
  [Python Setup and Usage]: https://docs.python.org/3/using
  [Python website]: https://www.python.org/
  [virtual environment]: https://realpython.com/what-is-pip

### Install with pip { data-toc-label="with pip" }

Zensical can be installed into a virtual environment with `pip`.

=== ":material-apple: macOS"
    Open up a terminal window and install Zensical by first setting up a virtual
    environment and then using `pip` to install the Zensical package into it:

    ``` sh
    python3 -m venv .venv
    source .venv/bin/activate
    pip install zensical
    ```

=== ":fontawesome-brands-windows: Windows"
    Open up a Command Window and install Zensical by first setting up a virtual
    environment and then using `pip` to install the Zensical package into it:

    ```
    python3 -m venv .venv
    .venv\Scripts\activate
    pip install zensical
    ```

=== ":material-linux: Linux"
    Open up a terminal window and install Zensical by first setting up a virtual
    environment and then using `pip` to install the Zensical package into it:

    ``` sh
    python3 -m venv .venv
    source .venv/bin/activate
    pip install zensical
    ```

  [upgrade to the next major version]: upgrade.md
  [Python Markdown]: https://python-markdown.github.io/
  [Pygments]: https://pygments.org/
  [Python Markdown Extensions]: https://facelessuser.github.io/pymdown-extensions/
  [Using Python's pip to Manage Your Projects' Dependencies]: https://realpython.com/what-is-pip/

### Install with uv { data-toc-label="with uv" }

If you are developing software using Python, chances are you're already using
[`uv`][uv] as a package manager, which has become popular in recent years. To
install Zensical with `uv`, use:

[uv]: https://docs.astral.sh/uv/

=== ":material-apple: macOS"

    ```
    uv init
    uv add zensical
    ```

=== ":fontawesome-brands-windows: Windows"

    ```
    uv init
    uv add zensical
    ```

=== ":material-linux: Linux"

    ```
    uv init
    uv add zensical
    ```


### Usage

The general command line syntax for Zensical is:

```sh
zensical COMMAND [OPTIONS] [ARGS]...
```

### Commands

- [`new`](new.md)
- [`build`](build.md)
- [`serve`](preview.md)

### Help

- General help: `zensical --help`
- Command-specific help: `zensical <command> --help`



## Preview

You can start a local web server to preview your documentation site as you write.
This allows you to view your site in a web browser without deploying it to a
remote server.

!!! warning "Use this command for preview only"

    Note that the web server that is built into Zensical is intended for preview
    purposes only. It is not designed for production deployment. We recommend to
    use a dedicated web server like nginx or Apache.

### Usage

```sh
zensical serve [OPTIONS]
```

This starts a local web server that serves your documentation site on
[localhost:8000][live preview]. As you make changes to source files, the browser
will automatically reload the page you're on.

  [live preview]: http://localhost:8000

### Options

The `serve` command accepts the following options:

| Option                     | Short | Description                                   |
| -------------------------- | ----- | --------------------------------------------- |
| --config-file              | -f    | Path to the config file to use.               |
| --open                     | -o    | Open preview in default browser               |
| --dev-addr &lt;IP:PORT&gt; | -a    | IP address and port (default: localhost:8000) |
| --help                     |       | Show a help message and exit.                 |


## Build

To build your documentation site, run `zensical build`.

[site_dir]: ../setup/basics.md#site_dir

### Usage

```sh
zensical build [OPTIONS] 
```

This will generate the static site in the configured [`site_dir`][site_dir],
with the default being `site`.

### Options

You can run `zensical build --help` to get command-line help for the `build`
command. It supports the following options:

| Option        | Short | Description                     |
| ------------- | ----- | ------------------------------- |
| --config-file | -f    | Path to the config file to use. |
| --clean       | -c    | Clean cache.                    |
| --help        |       | Show a help message and exit.   |