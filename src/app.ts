import express from 'express';
import { Request, Response } from './types';
import { NewFolderClass } from './new-folder';

const app = express();
const port = process.env.PORT || 3000;

// Middleware to parse JSON bodies
app.use(express.json());

// Initialize the new folder class
const newFolderInstance = new NewFolderClass();

// Example route
app.get('/', (req: Request, res: Response) => {
    res.send('Hello World!');
});

// Additional routes can be added here

app.listen(port, () => {
    console.log(`Server is running on http://localhost:${port}`);
});