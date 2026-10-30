---
tags:
  - contributing
  - documentation
  - user manual
  - markdown
---

# Markdown

Markdown is a lightweight markup language for authoring text such as technical
documentation, blog posts, etc. It was created in 2004 by John Gruber, with
help from Aaron Swartz. Compared to alternatives like HTML or XML, it is easy to
write and read, using a lean and intuitive syntax.

Zensical uses [Python Markdown] for compatibility with Material for MkDocs.
Python Markdown aims to be a faithful and extensible port of [John
Gruber's original Markdown syntax][gruber]. For the most part, it is fully compatible
with the original specification, supporting the core features like headings,
lists, links, blockquotes, and inline formatting (e.g., bold and italics).

Both Python Markdown itself and the [Python Markdown Extensions] that Zensical
also supports provide extensions to the core Markdown language to cater for the
needs of technical writers who want to produce clear, compelling, and visually
attractive documentation.

[Python Markdown]: https://python-markdown.github.io/
[Python Markdown Extensions]: https://facelessuser.github.io/pymdown-extensions/

## Learning Markdown

The [original description of Markdown][gruber] by John Gruber is still a good
starting point if you are unfamiliar with Markdown syntax. The [Markdown Guide]
is another great resource.

[Markdown Guide]: https://www.markdownguide.org/

## Page title

Zensical produces a page title for each page. It uses, in order of priority:

1. a title defined for the page in the `nav` configuration setting.
2. a title defined in the front-matter of the page.
3. a first-level Markdown header in the page content
4. the base name of the Markdown file

So, explicitly defined navigation takes the highest priority, followed by an
explicitly defined page title, followed by a title derived from the first level
one Markdown header. If none of these are available, the filename is used as a
fallback.

## Headers
```
# H1 Header
## H2 Header
### H3 Header
#### H4 Header
##### H5 Header
###### H6 Header
```

## Text formatting
```
**bold text**
*italic text*
***bold and italic***
~~strikethrough~~
`inline code`
```

## Links and images
```
[Link text](https://example.com)
[Link with title](https://example.com "Hover title")
![Alt text](image.jpg)
![Image with title](image.jpg "Image title")
```

## Lists
```
Unordered:
- Item 1
- Item 2
  - Nested item

Ordered:
1. First item
2. Second item
3. Third item
```

## Blockquotes
```
> This is a blockquote
> Multiple lines
>> Nested quote
```

## Code blocks
````
```javascript
function hello() {
  console.log("Hello, world!");
}
```
````

## Tables
```
| Header 1 | Header 2 | Header 3 |
|----------|----------|----------|
| Row 1    | Data     | Data     |
| Row 2    | Data     | Data     |
```

## Horizontal rule
```
---
or
***
or
___
```

## Task lists
```
- [x] Completed task
- [ ] Incomplete task
- [ ] Another task
```

## Escaping characters
```
Use backslash to escape: \* \_ \# \`
```

## Line breaks
```
End a line with two spaces  
to create a line break.

Or use a blank line for a new paragraph.
```